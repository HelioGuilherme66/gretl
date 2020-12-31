BINDIR = ../bin

TESTOUT = $(BINDIR)/gretl_test_scripts
MAKEOUT = $(BINDIR)/gretl_scripts_baseline
COMPARE = $(BINDIR)/gretl_output_compare
MEMTEST = $(BINDIR)/gretl_leak_check
PROFILE = $(BINDIR)/gretl_scripts_profile
REFACTOR = $(BINDIR)/refactor_data

.PHONY:

test:
	$(TESTOUT)
	$(REFACTOR)
	@cat fails
	@echo "=== `pwd` ===" > diffs
	$(COMPARE) >> diffs

new:
	$(MAKEOUT)

memtest:
	$(MEMTEST)

profile:
	$(PROFILE)

list:
	ls -1 *.inp > ps.list

clean:
	rm -rf ./newout
	rm -f diffs fails errs out
	rm -f *.log
	rm -f gpttmp* session.inp
