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

Change the scoring system by changing the final variables.

    public static final int MATCH = x;
    public static final int MISMATCH = x;
    public static final int INDEL = x;

Call the alignStrands method in main.

    alignStrands("GATTACA", "GCATGCU");

Compile and run the java file.

    javac NeedlemanWunsch.java
    java NeedlemanWunsch 

## License
This java implementation is released under the [MIT License](LICENSE).
