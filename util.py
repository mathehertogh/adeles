
DEBUG = False

def debug(msg):
	global DEBUG
	if DEBUG:
		print("    " + msg)