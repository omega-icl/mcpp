PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

all : test_MLP test_Torch test_MLPREG test_MLP_tangent test_MLP_adjoint test_MLP_grad test_MLP_copy test_MLP_ident test_MLP_home test_MLP_torch test_MLP_onnx test_MLPREG_torch test_MLP_reverse test_KAN test_FFKAN sim_NN

test_MLP : test_MLP.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_Torch : test_Torch.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_MLPREG : test_MLPREG.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_MLP_tangent : test_MLP_tangent.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_MLP_adjoint : test_MLP_adjoint.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_MLP_grad : test_MLP_grad.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_MLP_copy : test_MLP_copy.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_MLP_ident : test_MLP_ident.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_MLP_home : test_MLP_home.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_MLP_torch : test_MLP_torch.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_MLP_onnx : test_MLP_onnx.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_MLPREG_torch : test_MLPREG_torch.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_MLP_reverse : test_MLP_reverse.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_KAN : test_KAN.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_FFKAN : test_FFKAN.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

sim_NN : sim_NN.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

%.o: %.cpp
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $< -o $@

#####

graph : gendot gensvg

gendot :
	dot -Tpng -O $(shell ls *.dot)

gensvg :
	@for FILE in $(shell ls *.png); do \
	  echo $$FILE; \
	  convert -resize 70% $$FILE $$FILE; \
	  (display $$FILE &); \
	done

#####

clean :
	rm -f test_MLP test_Torch test_MLPREG test_MLP_tangent test_MLP_adjoint test_MLP_grad test_MLP_copy test_MLP_ident test_MLP_home test_MLP_torch test_MLP_onnx test_MLPREG_torch test_MLP_reverse test_KAN test_FFKAN sim_NN
	-(rm -f *.o *.out *.eps *.dot *.png *.lp *~)
