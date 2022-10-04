import os, importlib

handledphis = ["128", ""]

_filelist = [filename for filename in os.listdir("generated/") if ".py" in filename]

_primesizes = [name[6:-3] for name in _filelist if "primes" in name]

primesdict = {int(prime): getattr(importlib.import_module("generated." + filename[:-3]), "PRIMES" + prime) for filename in _filelist for prime in _primesizes if "primes" + prime in filename}

pmnsdicts = {int(prime + sphi): getattr(importlib.import_module("generated." + filename[:-3]), "pmns" + sphi + "dict") for filename in _filelist for prime in _primesizes for sphi in handledphis if "generated" + prime + "pmns" + sphi + ".py" in filename}

pmnsWBdicts = {int(prime + sphi): getattr(importlib.import_module("generated." + filename[:-3]), "pmns" + sphi + "dict") for filename in _filelist for prime in _primesizes for sphi in handledphis if "generatedwithbase" + prime + "pmns" + sphi + ".py" in filename}


handledphis += ["64"]
