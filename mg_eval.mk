# Extra target for building the comparator's MadGraph evaluator.
#
# Use this from a generated MadGraph subprocess directory together with the
# generated Makefile:
#
#   make -f Makefile -f /path/to/QCDHjjComparator/mg_eval.mk hgg_mg_eval
#
# The generated Makefile knows the MadLoop object list in $(PROCESS) and the
# loop/library flags in $(LINKLIBS).  This fragment only adds the missing
# executable target for hgg_mg_eval.f.

MG_EVAL ?= hgg_mg_eval

$(MG_EVAL): $(MG_EVAL).o $(PROCESS) makefile $(LIBS) | mg_eval_loop_cards
	$(FC) $(FFLAGS) -o $@ $(MG_EVAL).o $(PROCESS) $(LINKLIBS)

mg_eval_loop_cards:
	$(MAKE) -C $(ROOT)/Source treatCardsLoopNoInit

$(LIBDIR)libdhelas.$(libext):
	$(MAKE) -C $(ROOT)/Source libdhelas

$(LIBDIR)libmodel.$(libext):
	$(MAKE) -C $(ROOT)/Source libmodel
