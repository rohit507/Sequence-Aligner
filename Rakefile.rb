#!/usr/local/bin/rake

# This is the rakefile for our project, it berforms build
#  and testing automation in much the same way as a 
#  Makefile, but is significantly better suited for our 
#  chosen language and build structure.

##################  Constants  ###################

SC = "scalac"   # Scala Compiler
SE = "scala"    # Scala Executor

SF = "src"      # Source Folder
CF = "cls"      # Classfile Folder
LF = "lib"      # Libraries Folder
BF = "bin"      # Binaries Folder (where we keep symlinks to
                #                   useful applications) 

BioJ = "#{LF}/biojava.jar"              # Biojava libraries
BioAJ = "#{LF}/biojava-alignment.jar"   # Biojava Alignment libraries
ScaJ = "#{SF}/lib/scala-library.jar"    # Scala Libraries
AutJ = "#{LF}/autojar.jar"              # Autojar Libraries

CP = "#{BioJ}:#{BioAJ}:#{SF}:#{CF}"  # Class Path

MD = "#{SF}/Manifest.txt"   # Manifest Document

PS = "Proj4Scala.jar"   # Project 4 Scala Runtime Executable
PJ = "Proj4.jar"        # Project 4 Java Standalone Executable

Sources = FileList["#{SF}/*.scala"].gsub("#{SF}/","")
Classes = Sources.ext(".class")

SCOpts = "-deprecation -unchecked"  # Scala Compiler Options
JavaOpts = "-Xmx"                   # Java Runtime Options

# Task to take the large scala folder and regenerate the
#  entire useful directory structure. This is nice if you
#  wish to save space in transport yet still manage useful 
#  things with scala. 

task 'scala' => 'scala/' 

file 'scala/' => 'lib/scala.tgz' do
    if not Dir.exists?('scala') 
        sh "tar -zxvf lib/scala.tgz"
        sh "mv scala-2.9.2/ scala/"
    end
end

# our app with input of various sorts

task :run => "run:default"

namespace :run do

    task :default => :scala

    task :scala => PS do 
        sh "#{BF}/#{SE} " +
           " #{PS} amos/small/crp177.seq"
    end

end

# clean up to certain levels of detail, a spare clean for
#  all of the code related files, and a full clean that'll
#  expunge everything even vaguely dublicated so we can 
#  package things up sensibly small. 

task :clean => "clean:default"

namespace :clean do

    task :default => [:classes,:jars]

    task :full do
        rm_rf 'scala/'
    end

    task :classes do
        rm_rf 'cls/'
        mkdir 'cls'
    end

    task :jars do 
        rm_rf PS 
    end
end

# tasks to build various executables and bits of final output

task :build => "build:default"

namespace :build do

    task :default => :partial 

    task :partial => PS

end


file PS => [:scala,"#{CF}/Project4.class",BioJ] do
    cp BioJ , PS
    sh "jar umf #{MD} #{PS} -C #{CF} ."
end


# Generate all the rules for the class files
Sources.each do |sourceFile|
    
    s = sourceFile.gsub(".scala","")
    
    file "#{CF}/#{s}.class" =>
                ["#{SF}/#{s}.scala",BioJ] do
        sh "#{BF}/#{SC} #{SCOpts} -classpath #{CP} #{SF}/#{s}.scala -d #{CF}"
    end

end

# Use dependency rules for the sourcefiles to make sure
#  each set of class files is generated in the proper order
#  and we can assemble a coherent jar out of the thing. 

file "#{CF}/KmerTable.class" => ["scala"]


file "#{CF}/BioLibs.class" => ["scala",
                               "#{CF}/KmerTable.class"]

file "#{CF}/Project4.class" => ["scala",
                                "#{CF}/BioLibs.class",
                                "#{CF}/KmerTable.class"]

# FIXME TODO FIXME TODO FIXME TODO FIXME TODO
# Below this is the old rakefile commented out
#  none of this should be used by the time we're
#  done, it should all have been moved above this
#  bar or removed entirely.

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


