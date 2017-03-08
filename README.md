# Nucleotide Alignment - Needleman Wunsch
[![Code Climate](https://codeclimate.com/github/EscapingGithub/Needleman-Wunsch/badges/gpa.svg)](https://codeclimate.com/github/EscapingGithub/Needleman-Wunsch)
[![codebeat badge](https://codebeat.co/badges/676857fc-0080-4473-89af-2e2e09a34fcb)](https://codebeat.co/projects/github-com-escapinggithub-needleman-wunsch-master)

**Version:** 1.0.0

**Authors:** Andrew Quach and Tamir Enkhjargal

## Introduction

The repository contains a quick java implementation of aligning 
nucleotide sequences utilizing the Needleman-Wunsh algorithm. 
This was done as a project for Computer Science A as an example
of dynamic programming.

## Usage

Create a NeedlemanWunsch.java object. See constructors for more information.

    // Two strands
    NeedlemanWunsch n1 = new NeedlemanWunsch("ABC", "CAB");
    // Two strands + no mismatching allowed
    NeedlemanWunsch n2 = new NeedlemanWunsch("ABC", "CAB", false);
    // Two strands + different MATCH/MISMATCH/INDEL scoring system
    NeedlemanWunsch n3 = new NeedlemanWunsch("ABC", "CAB", 2, -2, -1);
    // Two strands + different MATCH/MISMATCH/INDEL scoring system + no mismatching allowed
    NeedlemanWunsch n4 = new NeedlemanWunsch("ABC", "CAB", 2, -2, -1, false);

Call methods on the object.

    // Prints out aligned strands and alignment score
    n1.printStrandInfo();
    // Get score
    n2.getScore();
    // Get solution matrix
    n3.getSolution();
    // Get alignedStrands;
    n4.getAlignedStrands();

Compile and run the java files. See the example Main.java.

    javac *.java
    java Main

## License
This java implementation is released under the [MIT License](LICENSE).
