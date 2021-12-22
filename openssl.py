#!usr/bin/python3.8

import subprocess
from subprocess import Popen, PIPE
import random
import base64

import json

class OpensslError(Exception):
  pass

def encrypt(plaintext, passphrase, cipher='aes-128-cbc'):

  pass_arg = 'pass:{0}'.format(passphrase)
  args = ['openssl', 'enc', '-' + cipher, '-base64', '-pass', pass_arg, '-pbkdf2']
	

  if isinstance(plaintext, str):
    plaintext = plaintext.encode('utf-8')
	

  result = subprocess.run(args, input=plaintext, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


  error_message = result.stderr.decode()
  if error_message != '':
    raise OpensslError(error_message)


  '''
	pipeline = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
	stdout, stderr = pipeline.communicate(plaintext)
	error_message = stderr.decode()
	if error_message != '':
			raise OpensslError(error_message)
	return stdout.decode()
	'''
  return result.stdout.decode()
	

def decrypt(ciphertext, passphrase, cipher='aes-128-cbc'):
    
  pass_arg = 'pass:{0}'.format(passphrase)
  args = ['openssl', 'enc', '-d', '-'+cipher, '-base64', '-pass', pass_arg, '-pbkdf2']

  if isinstance(ciphertext, str):
    ciphertext = ciphertext.encode('utf-8')

  result = subprocess.run(args, input=ciphertext, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  error_message = result.stderr.decode()
  if error_message != '':
    raise OpensslError(error_message)

  '''	
	pipeline = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
	stdout, stderr = pipeline.communicate(ciphertext)
	error_message = stderr.decode()
	if error_message != '':
			raise OpensslError(error_message)
	return stdout.decode()
	'''
  return result.stdout.decode()
  
def encode_base64(text):
	args = ['base64']
	pipeline = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
	if isinstance(text, str):
		text = text.encode('utf-8')
	stdout, stderr = pipeline.communicate(text)
	error_message = stderr.decode()
	if error_message != '':
		raise OpensslError(error_message)
	return stdout

def decode_base64(text):
	args = ['base64', '-d']
	pipeline = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
	if isinstance(text, str):
		text = text.encode('utf-8')
	stdout, stderr = pipeline.communicate(text)
	error_message = stderr.decode()
	if error_message != '':
		raise OpensslError(error_message)
	return stdout


def encrypt_pub(plaintext, pub):
	args = ['openssl', 'pkeyutl', '-encrypt', '-pubin', '-inkey', pub]
	if isinstance(plaintext, str):
		plaintext = plaintext.encode('utf-8')
	pipeline = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
	stdout, stderr = pipeline.communicate(plaintext)
	error_message = stderr.decode()
	if error_message != '':
		raise OpensslError(error_message)
	return encrypt_base64(stdout).decode()

def decrypt_pub(cryptedtext, privatekey='my_private_key'):
	tmp= decrypt_base64(cryptedtext)
	args = ['openssl','pkeyutl','-decrypt','-inkey',privatekey]
	pipeline2 = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
	stdout, stderr = pipeline2.communicate(tmp)
	print(stderr.decode())
	return stdout.decode()

def encrypt_hybrid(plaintext, pub):
	key = "a"+str(random.getrandbits(256))
	ret1 = encrypt_pub(key,pub)
	ret2 = encrypt(plaintext,key)
	return {"session_key":ret1,"payload":ret2}


def decrypt_hybrid(cryptedtext, privatekey='my_private_key', cipher='aes-128-cbc'):
    if(isinstance(cryptedtext,str)):
    	cryptedtext=json.loads(cryptedtext)
    key = decrypt_pub(cryptedtext["session_key"],privatekey)
    return decrypt(cryptedtext["payload"],key,cipher)


def sign(plaintext,privatekey='my_private_key'):
	args = ['openssl', 'dgst', '-sha256', '-sign', privatekey]
	if isinstance(plaintext, str):
		plaintext = plaintext.encode('utf-8')
	pipeline = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
	stdout, stderr = pipeline.communicate(plaintext)
	error_message = stderr.decode()
	if error_message != '':
		raise OpensslError(error_message)
	return encrypt_base64(stdout).decode()


def verify_sign(sign,plaintext,publickey='my_public_key'):
	sign = decrypt_base64(sign)
	tmp_file = open("verify_sign_tmp","wb")
	tmp_file.write(sign)
	tmp_file.close()
	args = ['openssl', 'dgst', '-sha256', '-verify', publickey, '-signature', "verify_sign_tmp"]
	if isinstance(plaintext, str):
		plaintext = plaintext.encode('utf-8')
	pipeline = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
	stdout, stderr = pipeline.communicate(plaintext)
	error_message = stderr.decode()
	if error_message != '':
		raise OpensslError(error_message)
	return stdout.decode()


def encrypt_with_H(plaintext, passphrase, cipher='aes-128-cbc'):

  pass_arg = 'pass:{0}'.format(passphrase)
  args = ['openssl', 'enc', '-' + cipher, '-md', 'sha256', '-pass', pass_arg, '-pbkdf2']
	

  if isinstance(plaintext, str):
    plaintext = plaintext.encode('utf-8')
	

  result = subprocess.run(args, input=plaintext, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  error_message = result.stderr.decode()
  if error_message != '':
    raise OpensslError(error_message)

  '''
	pipeline = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
	stdout, stderr = pipeline.communicate(plaintext)
	error_message = stderr.decode()
	if error_message != '':
			raise OpensslError(error_message)
	return stdout.decode()
	'''
  return result.stdout


def decrypt_with_H(ciphertext, passphrase, cipher='aes-128-cbc'):
    
  pass_arg = 'pass:{0}'.format(passphrase)
  args = ['openssl', 'enc', '-d', '-'+cipher, '-md', 'sha256', '-pass', pass_arg, '-pbkdf2']

  if isinstance(ciphertext, str):
    ciphertext = ciphertext.encode('utf-8')

  result = subprocess.run(args, input=ciphertext, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  error_message = result.stderr.decode()
  if error_message != '':
    raise OpensslError(error_message)

  '''	
	pipeline = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
	stdout, stderr = pipeline.communicate(ciphertext)
	error_message = stderr.decode()
	if error_message != '':
			raise OpensslError(error_message)
	return stdout.decode()
	'''
  return result.stdout

