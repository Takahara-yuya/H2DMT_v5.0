# 使用 Intel MPI C++ 编译器
CXX = mpiicpx

# 优化与并行编译选项
CXXFLAGS = -O3 -funroll-loops -march=native -flto -fopenmp

# 源码目录
SRC_DIR = src

# 目标可执行文件名
EXECUTABLE = H2DMT

# 只编译你指定的那几个源文件（请根据你的实际文件名调整，确保它们在 src/ 目录下）
SOURCES = \
    src/MTdata2D.cpp \
    src/MTfileread.cpp \
    src/MTfwd2D.cpp \
    src/MTinverse2D.cpp \
    src/MTmesh2D.cpp \
    src/H2DMT.cpp

# 生成对应的目标文件路径（在 src/obj/ 下）
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp, $(SRC_DIR)/obj/%.o, $(SOURCES))

# 如果你需要额外的头文件搜索路径，可以在这里用 -I 添加，例如：
# INCLUDES = -Iinclude

# 如果你有额外的链接库或链接选项，可以在这里添加，例如：
# LDFLAGS = -L/path/to/libs
# LDLIBS = -lmylib

all: $(EXECUTABLE)

# 链接目标文件，生成可执行文件
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

# 编译单个 .cpp 文件为 .o 文件，放入 src/obj/ 目录
$(SRC_DIR)/obj/%.o: $(SRC_DIR)/%.cpp | $(SRC_DIR)/obj
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDES)

# 如果 obj 目录不存在，则创建它
$(SRC_DIR)/obj:
	mkdir -p $(SRC_DIR)/obj

# 清理生成的 .o 文件和可执行文件
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
	rm -rf $(SRC_DIR)/obj