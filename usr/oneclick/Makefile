############################################################
# Copyright (C) 2001-2009, CompHEP Collaboration            
# Author: A. Sherstnev                                      
############################################################

.PHONY: all preparecomphepcode preparegencode relocate

all: preparecomphepcode preparegencode relocate

preparecomphepcode:
	make -C chepcode link
	@echo "CompHEP libraries have been created"

preparegencode:
	@if (test -f gencode/.done); then \
		make -C gencode; \
	else \
		make -C gencode gen; \
	fi
	@if (test ! -L models); then ln -s results/models .; fi
	@if (test ! -L tmp); then ln -s results/tmp .;	fi
	@echo "Numerical generator has been created"

relocate:
	make -C chepcode oneclicksetup

############################################################
.PHONY: clean totclean
clean:
	@make -C chepcode clean
	@make -C gencode clean
	@echo "both chepcode and gencode are clean"

totclean:
	@make -C chepcode totclean
	@make -C gencode totclean
	@make -C usercode clean
	@rm -f results/*.exe results/diag_view results/n_comphep results/.procent
	@rm -f session.dat .procent .*path
	@if (test -L models); then rm models; fi
	@if (test -L tmp); then rm tmp;	fi
	@echo "The oneclick code is fully clean"
