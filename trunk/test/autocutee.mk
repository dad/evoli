# cutee autogen: begin
CUTEE=./cutee
t_runners=t.protein.cutee.cc t.folder.cutee.cc t.population.cutee.cc 

%.cutee.cc: $(srcdir)/%.h
	$(CUTEE) -o $@ $<


runtest.cc: cutee
	$(CUTEE) -m -o $@

runtest-clean:
	rm -f autocutee.mk cutee *.o *.cutee.cc runtest runtest.cc
	touch autocutee.mk

EXTRA_PROGRAMS=runtest
runtest_SOURCES=runtest.cc $(test_files) $(t_runners)
runtest_DEPENDENCIES=cutee autocutee.mk
# cutee autogen: end 
