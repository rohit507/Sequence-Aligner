#!/usr/local/bin/rake

## This is the rakefile for our project, it berforms build
##  and testing automation in much the same way as a 
##  Makefile, but is significantly better suited for our 
##  chosen language and build structure.

## Will generate and run a app with a set of default perameters, allowing
##  you to overload one at a time
#task :run => :build do
#    Seq1 = ENV["Seq1"] || "public1-1.fa"
#    Seq2 = ENV["Seq2"] || "public1-2.fa"
#    match = ENV["Match"] || 1
#    mismatch = ENV["Mismatch"] || -1
#    gap = ENV["Gap"] || -1
#    puts "\n"
#    sh "java -jar Proj3.jar #{Seq1} #{Seq2} #{match} #{mismatch} #{gap}" 
#end

## tasks to create
#SC = "scalac"               # Scala Compiler  
#SE = "scala"                # Scala Executor

#ScF = "scala/bin/"          # Location Of Scala Binaries
#ScL = "scala/libcls/"       # Location Of Scala Library classes
#SF = "src"                  # Source Folder
#LF = "lib"                  # Library Folder
#CF = "cls"                  # Class Folder

#Lbj = "biojava.jar"         # BioJava Jar Name

#CP = "#{LF}/#{Lbj}:#{SF}:#{CF}" #ClassPath

## The various source files we use
#Sources = FileList['src/*.scala'].gsub("src/","").ext("")

#MF = "Manifest.txt"         # Manifest filename

#SJ = "Proj3Scala.jar"       # Scala Dependant Ouput
#PJ = "Proj3.jar"            # Project Output Jar File

## Run
## Clean
#task :clean do
#    rm_f SJ
#    rm_f PJ
#    rm_rf CF
#    mkdir CF
#end

## Build
#task :default => :build

## just a simple convinience alias
#task :build => PJ

## create the jar file that includes biojava and our program, and 
##  have our program set to run as default. 
#task SJ =>  ["#{LF}/#{Lbj}","#{SF}/#{MF}","#{CF}/Project3.class"] do
#    cp("#{LF}/#{Lbj}",SJ)
#    sh "jar umf #{SF}/#{MF} #{SJ} -C #{CF} ."
#end

#task PJ => SJ do
#    cp(SJ,PJ)
#    sh "jar uMf #{PJ} -C #{ScL} ."
#    sh "jar umf #{SF}/#{MF} #{PJ}"
#end

## FIXME : Add task to take care of unzipping the scala-library source files

## Make sure we've got proper dependancy chains for each source file
#
#task "#{SF}/Project3.scala" => "#{CF}/BioLib.class"
#task "#{SF}/BioLib.scala" => "#{CF}/Util.class"

## Generate a build function for each neccesary classfile
#Sources.each do |sourceFile|
#    rule "#{CF}/#{sourceFile}.class" =>
#                ["#{SF}/#{sourceFile}.scala","#{LF}/#{Lbj}"] do
#        sh "#{ScF + SC} -classpath #{CP} #{SF}/#{sourceFile}.scala -d #{CF}"
#    end
#end


