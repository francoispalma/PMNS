from generated.generated256pmns import pmnsdict as pmnsdict256
from generated.generated512pmns import pmnsdict as pmnsdict512
from generated.generated1024pmns import pmnsdict as pmnsdict1024
from generated.generated2048pmns import pmnsdict as pmnsdict2048
#from generated.generated4096pmns import pmnsdict as pmnsdict4096
#from generated.generated8192pmns import pmnsdict as pmnsdict8192
from generated.generated256pmns128 import pmns128dict as pmns128dict256
from generated.generated512pmns128 import pmns128dict as pmns128dict512
from generated.generated1024pmns128 import pmns128dict as pmns128dict1024
from generated.generated2048pmns128 import pmns128dict as pmns128dict2048
#from generated.generated4096pmns128 import pmns128dict as pmns128dict4096
#from generated.generated8192pmns128 import pmns128dict as pmns128dict8192
from generated.primes256 import PRIMES256
from generated.primes512 import PRIMES512
from generated.primes1024 import PRIMES1024
from generated.primes2048 import PRIMES2048
from generated.primes4096 import PRIMES4096
from generated.primes8192 import PRIMES8192

pmnsdict4096 = {}
pmns128dict4096 = {}
pmnsdict8192 = {}
pmns128dict8192 = {}

pmnsdicts = {256:pmnsdict256, 512:pmnsdict512, 1024:pmnsdict1024, 2048:pmnsdict2048, 4096:pmnsdict4096, 8192: pmnsdict8192,
        256128:pmns128dict256, 512128:pmns128dict512, 1024128: pmns128dict1024, 2048128:pmns128dict2048, 4096128: pmns128dict4096, 8192128:pmns128dict8192}

primesdict = {256: PRIMES256, 512: PRIMES512, 1024:PRIMES1024, 2048:PRIMES2048, 4096: PRIMES4096, 8192: PRIMES8192}

handledphis = ["128", "104", "64", ""]
