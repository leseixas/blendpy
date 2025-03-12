import os

LOGFILES = ['optimize.log', 'test_optimize.log', 'test_optimize_nostress.log']

def pytest_sessionfinish(session, exitstatus):
    print("\n\n  FINISHING TESTS")
    for logfile in LOGFILES:
        if os.path.exists(logfile):
            os.remove(logfile)
            print(f"    {logfile} has been deleted.")
