from generated.generated1024pmns import pmnsdict as pmnsdict1024
from generated.generated2048pmns import pmnsdict as pmnsdict2048
from generated.generated4096pmns import pmnsdict as pmnsdict4096
from generated.generated1024pmns128 import pmns128dict as pmns128dict1024
from generated.generated2048pmns128 import pmns128dict as pmns128dict2048
from generated.generated4096pmns128 import pmns128dict as pmns128dict4096
from generated.primes1024 import PRIMES1024
from generated.primes2048 import PRIMES2048
from generated.primes4096 import PRIMES4096

pmnsdicts = {1024:pmnsdict1024, 2048:pmnsdict2048, 4096:pmnsdict4096,
	1024128: pmns128dict1024, 2048128:pmns128dict2048, 4096128: pmns128dict4096}

primesdict = {1024:PRIMES1024, 2048:PRIMES2048, 4096: PRIMES4096}
