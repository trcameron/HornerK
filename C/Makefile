#####################################################
#   C Source Make File   							#
#####################################################
include ../make.inc
	
horner_test_acc: eft.h eft.c horner.h horner.c sum.h sum.c horner_test_acc.c
	$(CC) $(COPT) -o horner_test_acc eft.c horner.c sum.c horner_test_acc.c $(CFLAGS)
	
horner_test_time: eft.h eft.c horner.h horner.c sum.h sum.c horner_test_time.c
	$(CC) $(COPT) -o horner_test_time eft.c horner.c sum.c horner_test_time.c $(CFLAGS)
	
hornercmplx_test_acc: eft.h eft.c horner.h horner.c sum.h sum.c hornercmplx_test_acc.c
	$(CC) $(COPT) -o hornercmplx_test_acc eft.c horner.c sum.c hornercmplx_test_acc.c $(CFLAGS)
	
hornercmplx_test_time: eft.h eft.c horner.h horner.c sum.h sum.c hornercmplx_test_time.c
	$(CC) $(COPT) -o hornercmplx_test_time eft.c horner.c sum.c hornercmplx_test_time.c $(CFLAGS)
	
hornercomp_test: eft.h eft.c horner.h horner.c sum.h sum.c hornercomp_test.c
	$(CC) $(COPT) -o hornercomp_test eft.c horner.c sum.c hornercomp_test.c $(CFLAGS)
	
uninstall:
	@rm -f horner_test_acc
	@rm -f horner_test_time
	@rm -f hornercmplx_test_acc
	@rm -f hornercmplx_test_time
	@rm -f hornercomp_test