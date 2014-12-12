package Main;
import general.ExtendedReader;
import general.Functions;

import java.io.FileReader;
import java.util.Hashtable;
import Sequence.SamCoverage;
import Sequence.SequenceHandling;


public class Main {
	public static void main(String []args){
		int length = args.length;
		for (int i = 0; i < length; i++){
			args[i] = args[i].trim();
			System.out.print(args[i]+" ");
		}
		System.out.println();
		Hashtable<String,String> T = Functions.parseCommandLine(args);
		String program = Functions.getValue(T, "-p", "HELP").toUpperCase();
		if(program.indexOf("HELP") == 0){
			System.out.println("Following programs can be runned -p <option>");
//			System.out.println("csfasta");
//			CfastaSequences.run(T);
//			System.out.println("genome");
//			Genome.run(T);
//			System.out.println("databases");
//			Databases.run(T);
			System.out.println("sequencehandling");
			//SequenceHandling.run(T);
			System.out.println("sbatch");
			
		}
//
//		if(program.indexOf("CSFASTA") == 0){
//			CfastaSequences.run(T);
//		}

		if(program.indexOf("genomeCov".toUpperCase()) ==0 ){
			SamCoverage.findCoverage(T);
			return;
		}
		if(program.indexOf("SEQUENCEHANDLING".toUpperCase()) == 0){
			System.out.println("SequenceHandling");
			SequenceHandling.run(T);
		}
		if(program.indexOf("tableParsing".toUpperCase()) == 0){
			System.out.println("SequenceHandling");
			SequenceHandling.run(T);
		}


		if(program.indexOf("parseFiles".toUpperCase()) == 0){
			parseFiles();
		}
	}
	
	
	public static void parseFiles(){
		String dir = "/Users/johanreimegard/Vetenskap/Databases/acrocona/bowtie/stderr";
		String [] A = new String[6];
		A[0] = "20111027_11_15_45.stderr.txt";
		A[1] = "20111027_13_48_00.stderr.txt";
		A[2] = "20111027_13_53_48.stderr.txt";
		A[3] = "20111027_13_53_51.stderr.txt";
		A[4] = "20111027_13_53_52.stderr.txt";
		A[5] = "20111027_13_53_58.stderr.txt";
		System.out.print("RefFile\t");
		System.out.print("M\t");
		System.out.print("N\t");
		System.out.print("L\t");
		System.out.print("X\t");
		System.out.print("min\t");
		System.out.print("LP\t");
		System.out.print("MP\t");
		System.out.print("NFP\t");
		System.out.print("DSP\t");
		System.out.print("LR\t");
		System.out.print("MR\t");
		System.out.print("NFR\t");
		System.out.println("s");

		
		for(int i = 0; i < 30; i++){
			for(int j = 0; j < 6 ; j++){
				printInfoBowtieStdErr(dir+"/"+i+"_bowtie2_fastq_"+i+"_"+A[j]);
			}
		}
	}
	
	private static void printInfoBowtieStdErr(String file){
		int M, N,L,X,min;
		String	RefFile;
		if(file.compareTo("/Users/johanreimegard/Vetenskap/Databases/acrocona/bowtie/stderr/20_bowtie2_fastq_20_20111027_11_15_45.stderr.txt") != 0){
		try{
			
			ExtendedReader ER = new ExtendedReader(new FileReader(file));
			while(ER.more()){
//				bowtie2 --local -M 8 -N 2 -L 23 -X 1000  -t --score-min L,0,0.35 -x /bubo/proj/b2011098/private/bowtie2Ref/generated_0 -1 blend/blend.1.fastq -2 blend/blend.2.fastq -S /bubo/proj/b2011098/private/bowtie2/blend.1.fastq.blend.2.fastq.generated_0.28.sam
				String Line = ER.readLine();
				String[] temp = Line.split(" ");
				M = Integer.parseInt(temp[3]);
				N = Integer.parseInt(temp[5]);
				L = Integer.parseInt(temp[7]);
				X = Integer.parseInt(temp[9]);
				String[] t2 = temp[13].split("\\.");
				min = Integer.parseInt(t2[1]);
				String[] Refs = temp[15].split("/");
				RefFile= Refs[Refs.length-1];
//				Time loading reference: 00:00:00
				ER.skipLine();
//				Time loading forward index: 00:00:02
				ER.skipLine();
//				Time loading mirror index: 00:00:02
				if(N > 0)
					ER.skipLine();
//				Multiseed full-index search: 00:00:23
				ER.skipLine();
//				54000 reads; of these:
				ER.skipLine();
//				  54000 (100.00%) were paired; of these:
				ER.skipLine();
//				    53748 (99.53%) aligned concordantly 0 times
				ER.skipLine();
//				    245 (0.45%) aligned concordantly exactly 1 time
				int LP= ER.readInt();
				ER.skipLine();
//				    7 (0.01%) aligned concordantly >1 and <=8 times
				int MP = ER.readInt();
				ER.skipLine();
//				    0 (0.00%) aligned concordantly >8 times
				int NFP = ER.readInt();
				ER.skipLine();
//				    ----
				ER.skipLine();
//				    53748 pairs aligned concordantly 0 times; of these:
				ER.skipLine();
//				      9 (0.02%) aligned discordantly 1 time
				int DSP = ER.readInt();
				ER.skipLine();
//				    ----
				ER.skipLine();
//				    53739 pairs aligned 0 times concordantly or discordantly; of these:
				ER.skipLine();
//				      107478 mates make up the pairs; of these:
				ER.skipLine();
//				        106421 (99.02%) aligned 0 times
				ER.skipLine();
//				        1048 (0.98%) aligned exactly 1 time
				int LR = ER.readInt();
				ER.skipLine();
//				        9 (0.01%) aligned >1 and <=8 times
				int MR = ER.readInt();
				ER.skipLine();
//				0 (0.00%) aligned >8 times
				int NFR = ER.readInt();
				ER.skipLine();
//				1.46% overall alignment rate
				ER.skipLine();
//				Time searching: 00:00:27
				ER.skipLine();
//				Overall time: 00:00:27
				String t = ER.readLine();
				String[] time =t.split(":");
				int h = Integer.parseInt(time[1].trim());
				int m = Integer.parseInt(time[2].trim());
				int s = Integer.parseInt(time[3].trim());
				s = s + m*60+h*3600;
				
				System.out.print(RefFile+"\t");
				System.out.print(M+"\t");
				System.out.print(N+"\t");
				System.out.print(L+"\t");
				System.out.print(X+"\t");
				System.out.print(min+"\t");
				System.out.print(LP+"\t");
				System.out.print(MP+"\t");
				System.out.print(NFP+"\t");
				System.out.print(DSP+"\t");
				System.out.print(LR+"\t");
				System.out.print(MR+"\t");
				System.out.print(NFR+"\t");
				System.out.println(s);
			}

		}catch(Exception E){System.out.println(file);}
		
		}
		
		
	}


}

