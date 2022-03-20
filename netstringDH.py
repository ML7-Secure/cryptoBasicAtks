""" 
Python twisted helper v1.1. Copyright (C) LIP6

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; By running this program you implicitly agree
that it comes without even the implied warranty of MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE and you acknowledge that it may 
potentially COMPLETELY DESTROY YOUR COMPUTER (even if it is unlikely), 
INFECT IT WITH A VERY NASTY VIRUS or even RUN ARBITRARY CODE on it. 
See the GPL (GNU Public License) for more legal and technical details.
"""

from twisted.protocols.basic import NetstringReceiver


import hmac
from Crypto.Cipher import AES
from Crypto.Util import Counter

import random as rd
import subprocess
import json
from hashlib import sha256
from crypta import * ###

class ProtocolTransportMixin:
    """
    Act both as a protocol (w.r.t the outside world) AND a transport.
    """
    def __init__(self, protocolFactory, *args, **kwds):
        self.inner_protocol_factory = protocolFactory
        self.inner_protocol_args = args
        self.inner_protocol_kwds = kwds
        super().__init__()

        
        ## Encryption attributs ##
        self.Kaes_enc = bytes(16)
        self.Kiv_enc = 0 #int
        self.Kmac_enc = bytes(16)

        self.ctr_enc = Counter.new(128, initial_value=self.Kiv_enc)
        self.cipher_enc = AES.new(self.Kaes_enc, AES.MODE_CTR, counter=self.ctr_enc)
        
        ## Decryption attributs ##
        self.Kaes_dec = bytes(16)
        self.Kiv_dec = 0 #int
        self.Kmac_dec = bytes(16)
        
        self.ctr_dec = Counter.new(128, initial_value=self.Kiv_dec)
        self.cipher_dec = AES.new(self.Kaes_dec, AES.MODE_CTR, counter=self.ctr_dec)
        
        ## DH attributs ##
        self.dicoReceived = False
        self.dicoSent = False
        self.go_enc = False
        
        self.A = 0
        self.x = 0
        self.p = 0
        
        
        
    def _connect_inner_protocol(self):
        """
        Connects the inner protocol. Give it "self" as transport.
        """
        # build the inner protocol from the factory
        self.inner_protocol = self.inner_protocol_factory(*self.inner_protocol_args, **self.inner_protocol_kwds)

        # I'm the transport for the inner protocol
        inner_transport = self
        self.inner_protocol.makeConnection(inner_transport)

    def connectionMade(self):
        """
        I've been started up. Start the inner protocol.
        """
        super().connectionMade() # make connection in potential super-class
        self._connect_inner_protocol()

    # there is no dataReceived method: subclasses have to implement it

    def connectionLost(self, reason):
        """
        I've lost the connection to the outside world. Notify the inner protocol.
        """
        super().connectionLost(reason) # lose connection in potential super-class
        self.inner_protocol.connectionLost(reason) # lose connection in inner protocol


    #### now the transport part.
    def write(self, data):
        """
        Invoked by the inner protocol when it wants to send <data> to the outside world. 
        Intercept outgoing data, process, send to my own ("external") transport.
        """
        self.transport.write(data)

    def writeSequence(self, seq):
        """
        For CTRL+D
        """
        self.write(b''.join(seq))


    def loseConnection(self):
        """
        Inner protocol wants to abort. We abort.
        """
        self.transport.loseConnection()

    def getHost(self):
        """
        Inner protocol asks about the connected party. We forward the query...
        """
        return self.transport.getHost()

    def getPeer(self):
        return self.transport.getPeer()


class NetstringWrapperProtocol(ProtocolTransportMixin, NetstringReceiver):
    """
    This protocol receives netstrings. When a complete netstring is received,
    the payload is sent to the inner protocol. When the inner protocol sends
    bytes, they are wrapped inside a netstring and sent to the outside world.

    IMPLEMENTATION DETAILS :

    The connectionMade() method is inherited from ProtocolTransportMixin (both
    parent classes implement this method, but we inherit from
    ProtocolTransportMixin first, so it wins). It starts the inner proto AND
    call super().connectionMade(), which will resolve to
    NetstringReceiver.connectionMade()...

    The dataReceived() method is inherited from NetstringReceiver, because
    ProtocolTransportMixin does not implement it. It will call the
    stringReceived() method of this class once a netstring is received.
    """

    def aes_ctr_enc(self, data):
        """
        AES-128-CTR encryption
        using Cryptodome lib
        """
        return self.cipher_enc.encrypt(data)
    
        
    def aes_ctr_dec(self, data):
        """
        AES-128-CTR decryption
        using Cryptodome lib
        """
        return self.cipher_dec.decrypt(data)


    ################## Diffie-Hellman Methods ##################
    
    def dh_send(self): 
        """
        Initiate DH key exchange
        """
        # From the mussh pb_key
        self.p = 17125458317614137930196041979257577826408832324037508573393292981642667139747621778802438775238728592968344613589379932348475613503476932163166973813218698343816463289144185362912602522540494983090531497232965829536524507269848825658311420299335922295709743267508322525966773950394919257576842038771632742044142471053509850123605883815857162666917775193496157372656195558305727009891276006514000409365877218171388319923896309377791762590614311849642961380224851940460421710449368927252974870395873936387909672274883295377481008150475878590270591798350563488168080923804611822387520198054002990623911454389104774092183

        g = 8041367327046189302693984665026706374844608289874374425728797669509435881459140662650215832833471328470334064628508692231999401840332046192569287351991689963279656892562484773278584208040987631569628520464069532361274047374444344996651832979378318849943741662110395995778429270819222431610927356005913836932462099770076239554042855287138026806960470277326229482818003962004453764400995790974042663675692120758726145869061236443893509136147942414445551848162391468541444355707785697825741856849161233887307017428371823608125699892904960841221593344499088996021883972185241854777608212592397013510086894908468466292313

        q = 63762351364972653564641699529205510489263266834182771617563631363277932854227

        
        self.x = rd.randint(1,q)
        self.A = pow(g, self.x, self.p)

        with open('../login/A.txt', 'w') as f:
             f.write(str(self.A))

        #sign A
        args = ['openssl', 'dgst', '-hex', '-sha256', '-sign', '../login/private-key.pem', '../login/A.txt']
        result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        tmp = result.stdout.decode()
        signature = tmp.lstrip("EC-SHA256(../login/A.txt)= ")
        signature = signature.rstrip('\n').upper()
        #print('my signature : ', signature)

        dico = {'username': 'name.surname', 'A': self.A, 'signature': signature}
        myLog = json.dumps(dico).encode()
        #print('myLog : ', myLog)
        
        self.sendString(myLog)
        
    def dh_verif(self, dico):
        """
        Verify signature
        """

        B = dico['B']
        #print('B :', B)
        
        receivedSignature = dico['signature']
        #print('receivedSignature :', receivedSignature)
        
        S = str(self.A)+','+str(B)+','+"name.surname" #"<A>,<B>,<username>" 
        
        with open('../S.txt', 'w') as f:
             f.write(str(S))
        
        with open('../signature.bin', 'w') as f:
             f.write(receivedSignature)
             
        #signature from hex to bin (for openssl)
        args = ['xxd', '-r', '-p', '../signature.bin']
        result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        tmp = result.stdout
        #print('out xdd :',tmp)
             
        with open('../signature.bin', 'wb') as f:
             f.write(tmp)
        
        #verify signature
        args1 = ['openssl', 'dgst', '-sha256', '-verify', '../mussh_pkey.txt', '-signature', '../signature.bin', '../S.txt']
        result1 = subprocess.run(args1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        testSign = result1.stdout.decode()

        if "OK" in testSign:
            #print('\nVerification done :')
            print(testSign)
            K = pow(B, self.x, self.p)
            #print('K :', K)
            
            ## Decryption attributs S --> C ##
            
            a = (str(K)+'A').encode()
            self.Kaes_dec = sha256(a).digest()[:16] # [:16] to keep the 16 first bytes (bytes here)
            #print('self.Kaes_dec :', self.Kaes_dec)
            
            b = (str(K)+'B').encode()
            self.Kiv_dec = int(sha256(b).hexdigest()[:32], 16) # [:32] to keep the 16 first bytes (hex here)
            #print('self.Kiv_dec :', self.Kiv_dec)
            
            c = (str(K)+'C').encode()
            self.Kmac_dec = sha256(c).digest()[:16]
            #print('self.Kmac_dec :', self.Kmac_dec)
            
            self.ctr_dec = Counter.new(128, initial_value=self.Kiv_dec)
            self.cipher_dec = AES.new(self.Kaes_dec, AES.MODE_CTR, counter=self.ctr_dec)
            
            
            ## Encryption attributs C --> S ##
            
            d = (str(K)+'D').encode()
            self.Kaes_enc = sha256(d).digest()[:16] 
            #print('self.Kaes_enc :', self.Kaes_enc)
            
            e = (str(K)+'E').encode()
            self.Kiv_enc = int(sha256(e).hexdigest()[:32], 16)
            #print('self.Kiv_enc :', self.Kiv_enc)
            
            f = (str(K)+'F').encode()
            self.Kmac_enc = sha256(f).digest()[:16]
            #print('self.Kmac_enc :', self.Kmac_enc)

            self.ctr_enc = Counter.new(128, initial_value=self.Kiv_enc)
            self.cipher_enc = AES.new(self.Kaes_enc, AES.MODE_CTR, counter=self.ctr_enc)
            
            self.go_enc = True # First messages are not crypted
            self.communication()
            
        else:
            #print('\nVerification done :')
            print(testSign)
            exit('******************* '+testSign)
    
    ################## Diffie-Hellman Methods ##################
    
    def communication(self):
        self.write(b"test\n")
        
        # FEISTEL
        #self.feiselAttack()
        
    def feiselAttack(self):
        
        self.write(b"sus\n")
        
        #self.write(b"Single\n")
        self.write(b"Batch\n") 

        self.write(b"Enc\n")
        #self.write(b"Decryption\n")
        
        '''
        INPUTS (hex, 128 bits, 1 per line, 1024 lines max, empty line to end) :  
           0 >>> 7b2faf3a4787ea414434c1b211f5419e
           1 >>>  
        '''
               
    
    def stringReceived(self, data):
        """
        netstring received and decoded. Forward to inner proto.
        """
        #print('dataReceived : ', data)
        
        ## Catch DH parameters sent from the server ##
        if not self.dicoReceived:
           #self.inner_protocol.dataReceived(data)
           
           dico = json.loads(data)
           #print('dico B :', dico)
           self.dicoReceived = True
           
           self.dh_verif(dico)

        else:            
            ## DECRYPTION ##
            plain_data = self.aes_ctr_dec(data[:-32])
            #print('plain_data : ', plain_data)
            
            tag = hmac.new(self.Kmac_dec, plain_data, 'sha256').digest()
            
            # verify tag, bad tag => message ignored
            if tag != data[-32:] :
                return
                
            else:
                #print('plain_data : ', plain_data)
                self.inner_protocol.dataReceived(plain_data)

        #self.inner_protocol.dataReceived(data)
        
    def write(self, data):
        """
        Intercept outgoing data from inner protocol, wrap in netstring and
        send to external transport.
        """
        #print('SentData : ', data)
        
        ## Send DH parameters to the server ##
        if not self.dicoSent:
           self.dicoSent = True
           self.dh_send()
        
        #'''
        if self.go_enc:
            ## ENCRYPTION ##
            cipher_data = self.aes_ctr_enc(data)
            #print('cipher_data : ', cipher_data)
            tag = hmac.new(self.Kmac_enc, data, 'sha256').digest()
            auth_data = cipher_data+tag
            
            self.sendString(auth_data) # inherited from NetstringReceiver
        #'''
 
        #self.sendString(data) # inherited from NetstringReceiver
    
    

