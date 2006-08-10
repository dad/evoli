all: snpmisdg
#density
#tolerated-mutations
#test-fitness-density
#time-evol sequence-generator genebank-analyzer sequence-analyzer hillclimber misfold dgdist stabdist evolstab

run: snpmisdg
	./snp-mistrans-stability tr 10 6 0.0114 59.0 104.5 599 -5 -100 0 ../analysis/trs599ca6-snps.txt
#	./snp-fitness tr 1000 6 0.0114 59.0 104.5 599 -5 -100 0 ../analysis/trs599ca6-snps.txt
#	evolved-ddg-dist tr 1000 6 0.0114 59.0 104.5 599 -5 -100 0 ../analysis/trs599ca6-genes.txt
#	./misfold tr 1000 6 0.0114 59.0 104.5 599 -5 -100 0 ../analysis/trs599ca6-genes.txt 10
#	comparison tr 6 0 0.5 -5 -100 10000 ../seqs/start_sequence
#	density tr 6 0 0.002 -5 -100 ../seqs/start_sequence 
#	tolerated-mutations
#	test-fitness-density ../analysis/traa0ca6-genes.txt
#	misfold tr 1000 6 0 0.002 -5 -100 1000 ../analysis/traa0ca6-genes.txt 

tests: test-genetic-code test-genotype test-genotype-util \
	test-protein-folder test-fitness-evaluator \
	test-population test-translator test

#CPPFLAGS = -g -pedantic -Wall
#CPPFLAGS = -O3 -pedantic -Wall -pg
#CPPFLAGS = -O3 -pedantic -Wall
CPPFLAGS = -pedantic -Wall -fomit-frame-pointer -ffast-math -march=i686 -O3

#EXE = .exe
EXE = 

snpmisdg-obj = genetic-code.o protein-folder.o snp-mistrans-stability.o \
	translator.o fitness-evaluator.o

snpmisdg: $(snpmisdg-obj)
	g++ -o snp-mistrans-stability $(snpmisdg-obj)

snpstab_OBJ = genetic-code.o protein-folder.o snp-stability.o \
	translator.o

snpstab: $(snpstab_OBJ)
	g++ -o snp-stability $(snpstab_OBJ)

snpfit-obj = genetic-code.o protein-folder.o snp-fitness.o \
	translator.o fitness-evaluator.o

snpfit: $(snpfit-obj)
	g++ -o snp-fitness $(snpfit-obj)

evolddgdist_OBJ = genetic-code.o fitness-evaluator.o protein-folder.o evolved-ddg-dist.o \
	translator.o

evolddgdist: $(evolddgdist_OBJ)
	g++ $(CPPFLAGS) -o evolved-ddg-dist $(evolddgdist_OBJ)

tr-driver-nostr_OBJ = genetic-code.o fitness-evaluator.o translator.o \
	population.o translation-experiment-nostr.o tr-driver.o protein-folder.o  \

tr-driver-nostr: $(tr-driver-nostr_OBJ)
	g++ -o tr-driver-nostr $(tr-driver-nostr_OBJ)

tr-driver_OBJ = genetic-code.o fitness-evaluator.o translator.o \
	population.o translation-experiment.o tr-driver.o protein-folder.o  \

tr-driver: $(tr-driver_OBJ)
	g++ -o tr-driver $(tr-driver_OBJ)

comparison_OBJ = comparison.o genetic-code.o fitness-evaluator.o \
	population.o protein-folder.o  \
	translator.o

comparison: $(comparison_OBJ)
	g++ -o comparison $(comparison_OBJ)

neutral-ffold_OBJ = genetic-code.o neutral-ffold.o \
	protein-folder.o translator.o fitness-evaluator.o

neutral-ffold: $(neutral-ffold_OBJ)
	g++ -o neutral-ffold $(neutral-ffold_OBJ)

density_OBJ = genetic-code.o density.o population.o \
	protein-folder.o translator.o fitness-evaluator.o

density: $(density_OBJ)
	g++ -o density $(density_OBJ)

random-distributions_OBJ = genetic-code.o random-distributions.o population.o \
	protein-folder.o translator.o fitness-evaluator.o

random-distributions: $(random-distributions_OBJ)
	g++ -o random-distributions $(random-distributions_OBJ)

tolerated-mutations_OBJ = tolerated-mutations.o genetic-code.o \
	protein-folder.o translator.o

tolerated-mutations: $(tolerated-mutations_OBJ)
	g++ -o tolerated-mutations $(tolerated-mutations_OBJ)

fitness-advantage_OBJ = genetic-code.o fitness-advantage.o \
	protein-folder.o translator.o fitness-evaluator.o

fitness-advantage: $(fitness-advantage_OBJ)
	g++ -o fitness-advantage $(fitness-advantage_OBJ)

res-analyzer_OBJ = genetic-code.o res-analyzer.o \
	protein-folder.o translator.o

res-analyzer: $(res-analyzer_OBJ)
	g++ -o res-analyzer $(res-analyzer_OBJ)

evolstab_OBJ = genetic-code.o fitness-evaluator.o protein-folder.o evolstab.o \
	translator.o

evolstab: $(evolstab_OBJ)
	g++ $(CPPFLAGS) -o evolstab $(evolstab_OBJ)

stabdist_OBJ = genetic-code.o fitness-evaluator.o protein-folder.o stabdist.o \
	translator.o

stabdist: $(stabdist_OBJ)
	g++ $(CPPFLAGS) -o stabdist $(stabdist_OBJ)

dgdist_OBJ = genetic-code.o fitness-evaluator.o protein-folder.o stabdist.o \
	translator.o

dgdist: $(dgdist_OBJ)
	g++ $(CPPFLAGS) -o dgdist $(dgdist_OBJ)

misfold_OBJ = genetic-code.o fitness-evaluator.o \
	population.o misfold.o protein-folder.o  \
	translator.o

misfold: $(misfold_OBJ)
	g++ $(CPPFLAGS) -o misfold $(misfold_OBJ)

sequence-generator_OBJ = genetic-code.o sequence-generator.o \
	protein-folder.o translator.o

sequence-generator: $(sequence-generator_OBJ)
	g++ -o sequence-generator $(sequence-generator_OBJ)

sequence-analyzer_OBJ = genetic-code.o sequence-analyzer.o \
	protein-folder.o translator.o

sequence-analyzer: $(sequence-analyzer_OBJ)
	g++ -o sequence-analyzer $(sequence-analyzer_OBJ)

genebank-analyzer_OBJ = genetic-code.o genebank-analyzer.o \
	protein-folder.o translator.o fitness-evaluator.o

genebank-analyzer: $(genebank-analyzer_OBJ)
	g++ -o genebank-analyzer $(genebank-analyzer_OBJ)

hillclimber_OBJ = genetic-code.o hillclimber.o \
	protein-folder.o translator.o fitness-evaluator.o

hillclimber: $(hillclimber_OBJ)
	g++ -o hillclimber $(hillclimber_OBJ)

structure-printer_OBJ = structure-printer.o protein-folder.o \
	genetic-code.o translator.o

structure-printer: $(structure-printer_OBJ)
	g++ -o structure-printer $(structure-printer_OBJ)

time-evol_OBJ = genetic-code.o fitness-evaluator.o \
	population.o time-evol.o protein-folder.o \
        translator.o

time-evol: $(time-evol_OBJ)
	g++ -o time-evol $(time-evol_OBJ)

#
# All the test code
#

test_OBJ = test.o genetic-code.o \
	protein-folder.o translator.o fitness-evaluator.o

test: $(test_OBJ)
	g++ -o test $(test_OBJ)

test-fitness-density_OBJ = test-fitness-density.o genetic-code.o \
	protein-folder.o translator.o fitness-evaluator.o

test-fitness-density: $(test-fitness-density_OBJ) 
	g++ -o test-fitness-density $(test-fitness-density_OBJ)

test-translator_OBJ = translator.o genetic-code.o test-translator.o

test-translator: $(test-translator_OBJ)
	g++ -o test-translator $(test-translator_OBJ)

test-genetic-code_OBJ = genetic-code.o test-genetic-code.o

test-genetic-code: $(test-genetic-code_OBJ)
	g++ -o test-genetic-code $(test-genetic-code_OBJ)

test-genotype_OBJ = genetic-code.o test-genotype.o

test-genotype: $(test-genotype_OBJ)
	g++ -o test-genotype $(test-genotype_OBJ)

test-genotype-util_OBJ = genetic-code.o test-genotype-util.o \
	protein-folder.o translator.o fitness-evaluator.o

test-genotype-util: $(test-genotype-util_OBJ)
	g++ -o test-genotype-util $(test-genotype-util_OBJ)        

test-protein-folder_OBJ = genetic-code.o protein-folder.o \
	test-protein-folder.o

test-protein-folder: $(test-protein-folder_OBJ)
	g++ -o test-protein-folder $(test-protein-folder_OBJ)

test-fitness-evaluator_OBJ = genetic-code.o protein-folder.o translator.o \
	fitness-evaluator.o test-fitness-evaluator.o

test-fitness-evaluator: $(test-fitness-evaluator_OBJ)
	g++ -o test-fitness-evaluator $(test-fitness-evaluator_OBJ)

test-population_OBJ = genetic-code.o protein-folder.o fitness-evaluator.o \
	translator.o population.o test-population.o

test-population: $(test-population_OBJ)
	g++ -o test-population $(test-population_OBJ)

clean:
	rm -f $(transl-robustness_OBJ) transl-robustness$(EXE)             \
		$(misfold_OBJ) misfold$(EXE)									\
		$(sequence-generator_OBJ) sequence-generator$(EXE)           \
		$(genebank-analyzer_OBJ) genebank-analyzer$(EXE)            \
		$(hillclimber_OBJ) hillclimber-analyzer$(EXE)            \
		$(time-evol_OBJ) time-evol$(EXE)             \
		$(test-genetic-code_OBJ) test-genetic-code$(EXE)           \
		$(test-genotype_OBJ) test-genotype$(EXE)                   \
		$(test-protein-folder_OBJ) test-protein-folder$(EXE)       \
		$(test-fitness-evaluator_OBJ) test-fitness-evaluator$(EXE) \
		$(test-population_OBJ) test-population$(EXE)		     \
		$(test-translator_OBJ) test-translator$(EXE)                \
		$(test-genotype-util_OBJ) test-genotype-util$(EXE)        \
		$(dgdist_OBJ) dgdist$(EXE) 								\
		$(test_OBJ) test$(EXE) \
		$(transl-constraint_OBJ) transl-constraint$(EXE) \
		$(transl-acc_OBJ) transl-acc$(EXE) \
		$(res-analyzer_OBJ) res-analyzer$(EXE) \
		$(sequence-analyzer_OBJ) sequence-analyzer$(EXE) \
		$(test-fitness-density_OBJ) test-fitness-density$(EXE) \
		$(structure-printer_OBJ) structure-printer$(EXE) \
		$(transl-robonly_OBJ) transl-robonly$(EXE) \
		$(tolerated-mutations_OBJ) tolerated-mutations$(EXE) \
		$(density_OBJ) density$(EXE) \
		$(neutral-ffold_OBJ) neutral-ffold$(EXE) \
		$(comparison_OBJ) comparison$(EXE) \
		$(tr-driver_OBJ) tr-driver$(EXE) \
		$(tr-driver-nostr_OBJ) tr-driver-nostr$(EXE) \
		random-distributions$(EXE) $(random-distributions_OBJ) \
		evolved-ddg-dist$(EXE) $(evolddgdist_OBJ) \
		snp-stability$(EXE) $(snpstab_OBJ) \
		snp-mistrans-stability$(EXE) $(snpmisdg-obj) \
		snp-fitness$(EXE) $(snpfit-obj)

