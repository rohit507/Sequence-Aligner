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
TF = "tmp"      # Temp Folder
BF = "bin"      # Binaries Folder (where we keep symlinks to
                #                   useful applications) 

BioJ = "#{LF}/biojava.jar"              # Biojava libraries
BioAJ = "#{LF}/biojava-alignment.jar"   # Biojava Alignment libraries
ScaJ = "#{SF}/lib/scala-library.jar"    # Scala Libraries
AutJ = "#{LF}/autojar.jar"              # Autojar Libraries
ColtJ = "#{LF}/colt.jar"                # Colt Libraries
ColtCj = "#{LF}/colt-concurrent.jar"    # Concurrent Colt Libs
TroveJ = "#{LF}/trove.jar"              # Trobe Libraries

CP = "#{TroveJ}:#{SF}:#{CF}"  # Class Path

MD = "#{SF}/Manifest.txt"   # Manifest Document

PS = "Proj4Scala.jar"   # Project 4 Scala Runtime Executable
PJ = "Proj4.jar"        # Project 4 Java Standalone Executable

Sources = FileList["#{SF}/*.scala"].gsub("#{SF}/","")
Classes = Sources.ext(".class")

SCOpts = "-deprecation -unchecked -optimize"  # Scala Compiler Options
JavaOpts = "-Xmx2G"                           # Java Runtime Options
CmdLineArgs = ""                              # Runtime COmmand Line Options

ENV['JAVA_OPTS'] = JavaOpts

SmallInput = "amos/small/crp177.seq"    # The Short Test Sequence
LargeInput = "amos/c_ruddii.seq"        # The Long Test Sequence

# Utility function for time displays

def print_time_diff(diff)
    hrs = (diff / (60 * 60)).floor
    diff = diff % (60 * 60)
    min = (diff / 60).floor
    diff = diff % 60
    sec = diff.floor
    milli = ((diff % 1) * 1000).floor
    return "#{hrs}h:#{min}m:#{sec}s:#{milli}ms"
end

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
        act = ENV['args']||"" 
        sh "#{BF}/#{SE} #{PS} -i #{SmallInput} #{CmdLineArgs} #{act} "
    end

    task :large => PS do
        act = ENV['args']||"" 
        sh "#{BF}/#{SE} #{PS} #{act} -i #{LargeInput} #{CmdLineArgs} #{act}"
    end
end

# Run the full pipeline in a sensible fashion

task :pipeline => "pipeline:default"

namespace :pipeline do

    task :amos => "amos:default"

    namespace :amos do

        task :default do
            run_amos_pipe(SmallInput,true,true)
        end

        task :large do
            run_amos_pipe(LargeInput,true,false)
        end

        def run_amos_pipe(file,time,cat)
            rm_rf(TF)
            mkdir(TF)
            if File.exists?(file) then
                cp(file,TF)
            else 
                raise "No #{file} exists."
            end
            raw = file.pathmap("%n")
            seq = "#{TF}/" + raw.ext("seq")
            bnk = "#{TF}/" + raw.ext("bnk")
            fst = "#{TF}/" + raw.ext("fasta")
            strTime = Time.now()
            sh("#{BF}/toAmos_new -s #{seq} -b #{bnk}")
            bnkTime = Time.now()
            sh("#{BF}/hash-overlap #{bnk} -B -x 0.04 -o 40")
            ovrTime = Time.now()
            sh("#{BF}/tigger -b #{bnk}")
            tigTime = Time.now()
            sh("#{BF}/make-consensus -e 0.04 -o 40 -B -b #{bnk}")
            cnsTime = Time.now()
            sh("#{BF}/bank2fasta -b #{bnk} > #{fst}")
            fstTime = Time.now()

            if cat then
                sh("cat #{fst}")
            end

            if time then
                puts ""
                puts "============ Time Taken ============="
                puts "Total Time            : #{print_time_diff(fstTime - strTime)}"
                puts "  Bank Creation Time  : #{print_time_diff(bnkTime - strTime)}"
                puts "  Overlap Time        : #{print_time_diff(ovrTime - bnkTime)}"
                puts "  Contigger Time      : #{print_time_diff(tigTime - ovrTime)}"
                puts "  Consensus Time      : #{print_time_diff(cnsTime - tigTime)}"
                puts "  Fasta Creation Time : #{print_time_diff(fstTime - cnsTime)}"
                puts ""
            end
        end
    end

    task :project => "project:default"

    namespace :project do

        task :default => PS do
            run_project_pipe(SmallInput,true,true)
        end

        task :large => PS do
            run_project_pipe(LargeInput,true,false)
        end

        def run_project_pipe(file,time,cat)
            rm_rf(TF)
            mkdir(TF)
            if File.exists?(file) then
                cp(file,TF)
            else 
                raise "No #{file} exists."
            end
            raw = file.pathmap("%n")
            seq = "#{TF}/" + raw.ext("seq")
            ovl = "#{TF}/" + raw.ext("ovl")
            bnk = "#{TF}/" + raw.ext("bnk")
            fst = "#{TF}/" + raw.ext("fasta")
            strTime = Time.now()
            sh("#{BF}/toAmos_new -s #{seq} -b #{bnk}")
            bnkTime = Time.now()
                # Place your overlapper here
            sh("#{BF}/#{SE} #{PS} -i #{seq} -o #{ovl}" +
                " #{CmdLineArgs} #{ENV['args']}")
            ovrTime = Time.now()
            sh("#{BF}/bank-transact -b #{bnk} -m #{ovl}")
            trnTime = Time.now()
            sh("#{BF}/tigger -b #{bnk}")
            tigTime = Time.now()
            sh("#{BF}/make-consensus -e 0.04 -o 40 -B -b #{bnk}")
            cnsTime = Time.now()
            sh("#{BF}/bank2fasta -b #{bnk} > #{fst}")
            fstTime = Time.now()

            if cat then
                sh("cat #{fst}")
            end

            if time then
                puts ""
                puts "============ Time Taken ============="
                puts "Total Time               : #{print_time_diff(fstTime - strTime)}"
                puts "  Bank Creation Time     : #{print_time_diff(bnkTime - strTime)}"
                puts "  Overlap Time           : #{print_time_diff(ovrTime - bnkTime)}"
                puts "  Bank Transaction Time  : #{print_time_diff(trnTime - ovrTime)}"
                puts "  Contigger Time         : #{print_time_diff(tigTime - trnTime)}"
                puts "  Consensus Time         : #{print_time_diff(cnsTime - tigTime)}"
                puts "  Fasta Creation Time    : #{print_time_diff(fstTime - cnsTime)}"
                puts ""
            end
        end
    end

end

# clean up to certain levels of detail, a spare clean for
#  all of the code related files, and a full clean that'll
#  expunge everything even vaguely dublicated so we can 
#  package things up sensibly small. 

task :clean => "clean:default"

namespace :clean do

    task :default => [:classes,:jars,:temp]

    task :full do
        rm_rf 'scala/'
    end

    task :classes do
        rm_rf 'cls/'
        mkdir 'cls'
    end

    task :temp do
        rm_rf 'tmp/'
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
    cp TroveJ , PS
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

file "#{CF}/ObjectStore.class" => ["scala"]

file "#{CF}/KmerTable.class" => ["scala",
                                 "#{CF}/ObjectStore.class"]


file "#{CF}/BioLibs.class" => ["scala",
                               "#{CF}/KmerTable.class",
                               "#{CF}/ObjectStore.class"]

file "#{CF}/Project4.class" => ["scala",
                                "#{CF}/BioLibs.class",
                                "#{CF}/KmerTable.class",
                                "#{CF}/ObjectStore.class"]

