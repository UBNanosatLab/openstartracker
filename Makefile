.PHONY: all ost beast tests test clean production

TESTDIR ?= science_cam_may8_0.05sec_gain40

all: ost tests

ost:
	$(MAKE) -C ost

beast: ost

tests:
	$(MAKE) -C tests

test: tests
	cd tests && ./unit_test.sh -e $(TESTDIR)

production:
	./build_production.sh $(TESTDIR)

clean:
	$(MAKE) -C tests clean
	$(MAKE) -C beast clean
