SHELL := /bin/bash 
notebooks := $(wildcard *-*.ipynb)
notebooks_executed := $(foreach d,$(notebooks), notebooks4pdf/$(d) )

.PHONY: clean pdf notebooks dirs clean_output

clear_output: 
	@echo CLEANING...
	./clear_output.sh
	@echo OK


notebooks4pdf/%.ipynb : %.ipynb
	@echo CLEANING OUTPUT of $< 
	@jupyter nbconvert --ClearOutputPreprocessor.enabled=True --clear-output $<
	@echo executing notebook  $< and writing it to $@
	@time jupyter nbconvert  --to notebook --execute  $< --output $@  >& /tmp/nbconvert.log

pdf: $(notebooks_executed)  latex_template2.tplx 
	cd notebooks4pdf; python3 -m bookbook.latex --template ../latex_template2.tplx --pdf && mv combined.pdf ../dynamical_systems.pdf


clean:
	@rm -fv combined.*   *log *aux *tex
	@rm -frv notebooks4pdf/ html/

dirs:
	@test -d notebooks4pdf || mkdir -v notebooks4pdf && echo notebooks4pdf exists
	@test -L notebooks4pdf/images || ln -sv ../images notebooks4pdf/images && echo images link exists  

log1 := $(shell test -d notebooks4pdf || mkdir -v notebooks4pdf)
log2 := $(shell test -L notebooks4pdf/images || ln -sv ../images notebooks4pdf/images)
