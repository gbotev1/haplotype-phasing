# Define compiler variables
JC = javac
JVM = java
MAIN = Tester
JFLAGS = -Xlint:all -g -d bin
CP = -cp "/Users/gbotev/Documents/GitHub Repositories/haplotype-phasing/bin:/Users/gbotev/Documents/GitHub Repositories/haplotype-phasing/external_jars/guava-23.0.jar"
SOURCEPATH = src/*.java
OBJECTPATH = bin/*.class

# Default compilation
.SUFFIXES: .java .class

default: .java.class

.java.class:
	@$(JC) $(JFLAGS) $(CP) $(SOURCEPATH) 

clean:
	@$(RM) $(OBJECTPATH)

run:
	@$(JVM) $(CP) $(MAIN) $(args)