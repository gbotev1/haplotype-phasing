# Define compiler variables
JC = javac
JVM = java
MAIN = src.main.Tester
JFLAGS = -Xlint:all -g:none -d bin
CP = -cp "/Users/gbotev/Documents/GitHub Repositories/haplotype-phasing/bin:/Users/gbotev/Documents/GitHub Repositories/haplotype-phasing/external_jars/guava-23.0.jar"
SOURCEPATH = src/main/*.java
OBJECTPATH = bin/src/main/*.class

# Default compilation
.SUFFIXES: .java .class

default: .java.class

.java.class:
	@$(JC) $(JFLAGS) $(CP) $(SOURCEPATH) 

clean:
	@$(RM) $(OBJECTPATH)

run:
	@$(JVM) $(CP) $(MAIN)
