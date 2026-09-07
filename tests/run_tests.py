# zero-dependency test runner: imports tests/test_*.py and runs every test_ function
import importlib.util
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (ROOT, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)


def main():
    passed = failed = 0
    for fn in sorted(os.listdir(HERE)):
        if not (fn.startswith("test_") and fn.endswith(".py")):
            continue
        spec = importlib.util.spec_from_file_location(fn[:-3], os.path.join(HERE, fn))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        for name in sorted(dir(mod)):
            if not name.startswith("test_"):
                continue
            try:
                getattr(mod, name)()
                print(f"  PASS  {fn}::{name}")
                passed += 1
            except Exception:
                print(f"  FAIL  {fn}::{name}")
                traceback.print_exc()
                failed += 1
    print(f"\n{passed} passed, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
