# variable definition
CC		:= g++
SRC_DIR		:= src
OBJ_DIR		:= obj
EXEC_DIR	:= bin
SRC			:= $(wildcard $(SRC_DIR)/*.cc)
OBJ			:= $(addprefix $(OBJ_DIR)/, $(notdir $(SRC:.cc=.o)))
EXEC		:= $(EXEC_DIR)/rlacl
DEBUG		:= #-g3
OPT			:= -O3
WARN		:= #-Wall
MISC		:= -std=c++11 -c -fmessage-length=0
CFLAG		:= $(OPT) $(DEBUG) $(WARN) $(MISC)

INC			:= -I /opt/homebrew/Cellar/boost/1.79.0_2/include #/opt/local/include/libxml2
LIB			:= -lxml2

RM			:= rm -rf

# call all
all: dir $(OBJ)
	@echo "Building target: $(EXEC)"
	g++ -o $(EXEC) $(OBJ) $(LIB)
	@echo "Finished building target: $(EXEC)"
	@echo ""
	
dir:
	mkdir -p obj
	mkdir -p bin

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	@echo "Building file: $<"
	g++ $(INC) $(CFLAG)	-o "$@" "$<"
	@echo "Finished building file: $<"
	@echo
	
# call clean
clean:
	$(RM) $(OBJ) $(EXEC)
