import sys
from typing import Any
from fable_modules.fable_library.list import (of_array, FSharpList)
from fable_modules.fable_library.types import Array
from fable_modules.fable_library.util import (uncurry, to_enumerable)
from helpers import (initialize_context, run, dotnet, npm, run_parallel, run_or_default)

initialize_context()

client_path: str = Fake_IO_Path_getFullName("src/Client")

deploy_path: str = Fake_IO_Path_getFullName("deploy")

client_tests_path: str = Fake_IO_Path_getFullName("tests/Client")

def _arrow0(_arg: Any) -> None:
    Fake_IO_Shell_cleanDir(deploy_path)
    run(uncurry(2, dotnet), "fable clean --yes", client_path)


Fake_Core_TargetModule_create("Clean", _arrow0)

def _arrow1(_arg: Any) -> None:
    run(uncurry(2, npm), "install", ".")


Fake_Core_TargetModule_create("InstallClient", _arrow1)

def _arrow2(_arg: Any) -> None:
    run_parallel(to_enumerable([("client", dotnet("fable -o output -s --run npm run build")(client_path))]))


Fake_Core_TargetModule_create("Bundle", _arrow2)

def _arrow3(_arg: Any) -> None:
    run_parallel(to_enumerable([("client", dotnet("fable watch -o output -s --run npm run start")(client_path))]))


Fake_Core_TargetModule_create("Run", _arrow3)

def _arrow4(_arg: Any) -> None:
    run_parallel(to_enumerable([("client", dotnet("fable watch -o output -s --run npm run test:live")(client_tests_path))]))


Fake_Core_TargetModule_create("RunTests", _arrow4)

def _arrow5(_arg: Any) -> None:
    run(uncurry(2, dotnet), "fantomas . -r", "src")


Fake_Core_TargetModule_create("Format", _arrow5)

dependencies: FSharpList[str] = of_array([Fake_Core_TargetOperators_op_EqualsEqualsGreater(Fake_Core_TargetOperators_op_EqualsEqualsGreater("Clean", "InstallClient"), "Bundle"), Fake_Core_TargetOperators_op_EqualsEqualsGreater(Fake_Core_TargetOperators_op_EqualsEqualsGreater("Clean", "InstallClient"), "Run"), Fake_Core_TargetOperators_op_EqualsEqualsGreater("InstallClient", "RunTests")])

def main(args: Array[str]) -> int:
    return run_or_default(args)


if __name__ == "__main__":
    main(sys.argv[1:])


