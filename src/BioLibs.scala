/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/

import java.io._
import scala.collection._

object BioLibs {
    // Read in a fasta file dispatching each sequence as it's
    //  read through the input function. 
    def readSeq( file : String, act : (Int,String) => _) : Int = {
        var i = 1
        var s = ""
        var in = new BufferedReader(new FileReader(file))
        var line = in.readLine()

        if (!line.startsWith(">"))
            throw new Exception("Invalid Sequence File: " + file)

        line = in.readLine()
        while (line != null) {
            if (line.startsWith(">")) {
                act(i,s)
                i += 1
                s = ""
            } else {
                s += line
            }
            line = in.readLine()
        }
    
        act(i,s)

        return i
    }

    // Generates an integer unique to the first 16 bases of the
    //  of the input sequence, useful for a number of operations
    def seqHash( seq : String) : Int = {

        var h : Int = 0

        for (i <- 0 to 16) {
            var c = Character.toUpperCase(seq.charAt(i))

            h = h << 2 

            c match {
                case 'A' => h ^= 0
                case 'C' => h ^= 1
                case 'T' => h ^= 2
                case 'G' => h ^= 3
                case  _  => println("Char \'" + c + "\' is not " +
                                " acceptable input, Found in \"" +
                                seq + "\"")
            }
        }

        return h
    }

    // generates a mutable hashSet of strings 
    def generateKmerSet( k : Int )( s : String ) 
                        : (Int,Set[String]) = {
        var o = new mutable.HashSet[String]()

        for( i <- 0 to (s.length - k)) {
            o += s.substring(i,i+k)
        }
        return (k,o)
    }

    
}
