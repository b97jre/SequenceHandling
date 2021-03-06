package Sequence;

import general.Functions;
import general.GeneticCode;
import general.IOTools;
import general.RNAfunctions;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;

import general.ExtendedReader;
import general.ExtendedWriter;




public class FastaSequences extends ArrayList <FastaSequence> implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private String GenomeName;
	String Name;
	private boolean verbose;
	ArrayList <Pair> pairs;


	public FastaSequences(){
		verbose = false;
	}

	public FastaSequences(String file){
		this.verbose = false;
		readFastaFile(file);
	}

	public FastaSequences(String file,boolean verbose){
		this.verbose = verbose;
		readFastaFile(file);
	}

	public FastaSequences(String dir, String file){
		this.Name = file;
		verbose = false;
		parseAllFasta(file,dir);
	}

	FastaSequences(String dir, String file, int nrOfExperiments){
		verbose = false;
		readFastaFile(dir,file, nrOfExperiments);
	}


	public Hashtable<String, FastaSequence> convertToHashTable2(){
		Hashtable<String, FastaSequence> HT = new Hashtable<String, FastaSequence>();
		//System.out.println(this.size());
		for(int i = 0; i< this.size(); i++){
			String Name = this.get(i).Name;
			if(Name.indexOf(" ") >-1)
				Name = this.get(i).Name.split(" ")[0];
			HT.put(Name, this.get(i));
			//System.out.println(Name);
		}

		return HT;
	}



	public static FastaSequences getUniqueTranscripts(FastaSequences A, FastaSequences B,String DifferenceOut){
		FastaSequences C = new FastaSequences();
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(DifferenceOut));
			Hashtable<String, FastaSequence> HT = A.convertToHashTable2();

			int progress = 1;
			System.out.println(B.size());
			for(int i = 0; i < B.size();i++){
				String Name = B.get(i).Name;
				if(Name.indexOf(" ") >-1)
					Name = Name.split(" ")[0];			
				if(HT.containsKey(Name)){
					if(HT.get(Name).getLength() != B.get(i).getLength()){
						C.add(B.get(i));
						EW.println(Name+"\t"+HT.get(Name).getLength()+"\t"+B.get(i).getLength()+"\t"+(B.get(i).getLength() - HT.get(Name).getLength())+"\t"+((100*(B.get(i).getLength() - HT.get(Name).getLength()))/HT.get(Name).getLength()));
					}
				}else{
					C.add(B.get(i));
				}
				if(i> B.size()/100*progress){
					System.out.print(".");
					progress++;
				}
			}
			EW.flush();
			EW.close();
		}catch(Exception E){
			E.printStackTrace();
		}

		return C;
	}

	public static FastaSequences getExtension(FastaSequences A, FastaSequences B,String DifferenceOut){
		FastaSequences C = new FastaSequences();
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(DifferenceOut));
			Hashtable<String, FastaSequence> HT = A.convertToHashTable2();

			int progress = 1;
			System.out.println(B.size());
			for(int i = 0; i < B.size();i++){
				String Name = B.get(i).Name;
				if(Name.indexOf(" ") >-1)
					Name = Name.split(" ")[0];			
				if(HT.containsKey(Name)){
					if(HT.get(Name).getLength() != B.get(i).getLength()){
						int[] seq1 = HT.get(Name).getSequence();
						int[] seq2 = B.get(i).getSequence();
						int[] extension = RNAfunctions.getExtension(seq1,seq2);
						if(extension != null){
							EW.println(Name+"_"+RNAfunctions.getBeforeOrAfter(seq1,seq2));
							EW.println(RNAfunctions.DNAInt2String(extension));
						}else{
							System.out.println("Something strange with "+Name);
						}
					}
				}
				if(i> B.size()/100*progress){
					System.out.print(".");
					progress++;
				}
			}
			EW.flush();
			EW.close();
		}catch(Exception E){
			E.printStackTrace();
		}

		return C;
	}



	public static FastaSequences getMutualTranscripts(FastaSequences A, FastaSequences B){
		FastaSequences C = new FastaSequences();

		for(int i = 0; i < A.size();i++){
			int j = 0;
			while(j < B.size() && !B.get(j).sameSequence(A.get(i).Sequence)){
				j++;
			}
			if(j < B.size())
				C.add(B.get(j));
		}
		return C;
	}




	public int findDoubleMutations(int[][][] doubleMutations, int[] refSeq){
		int total = 0;
		for(int i = 0; i < this.size(); i++){
			if(this.get(i).hasTwoMutations(doubleMutations, refSeq)){
				total = total + 1 ;
			}
		}
		return total;

	}

	public  void printLength(){
		for(int i = 0; i < this.size(); i++){
			System.out.println(this.get(i).getLength());
		}
	}


	public void printGCcontent(){
		for(int i = 0; i < this.size(); i++){
			System.out.println(this.get(i).getName()+"\t"+this.get(i).getStringSequence()+"\t"+this.get(i).getGCcontent());
		}
	}


	public static void getORFs(String dir,String fastaFile){
		getORFs(dir+"/"+fastaFile);
	}



	public void printDistribution(String gmapperFile, int cutoff){
		String[] names = new String[this.size()];
		int []counts = new int[this.size()];
		int[] lengthDist = new int[35];

		for(int i = 0; i < this.size(); i++){
			System.out.println();
			int[][] dist = this.get(i).getDistr(gmapperFile);
			lengthDist = this.get(i).getSizeDistribution(lengthDist, dist);
			names[i] = this.get(i).getName();
			System.out.println(names[i]);
			System.out.println(this.get(i).getStringSequence());

			counts[i] = this.get(i).printDist(dist, cutoff);
			System.out.println();
		}


		System.out.println("Name\tcount");
		for(int i = 0; i < this.size(); i++){
			System.out.println(names[i]+"\t"+counts[i]);
		}


		System.out.println("size\tcount");
		for(int i = 0; i <	lengthDist.length; i++){
			System.out.println(i+"\t"+lengthDist[i]);
		}

	}

	int findMutants(String fileName){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			int count = 100000;
			while(ER.more()){

				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					FastaSequence a = getFasta(ER);
				}
				if(this.size() > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 100000;
				}
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
		return 1;
	}


	public int findSingleMutations(int[][] singleMutation, int[] refSeq){
		int total = 0;
		for(int i = 0; i < this.size(); i++){
			if(this.get(i).hasSingleMutation(singleMutation, refSeq)){
				total = total + 1 ;
			}
		}
		return total;

	}	


	private void addFasta(ExtendedReader ER, int nrOfExperiments){
		String InfoLine = ER.readLine();
		//>853_36_1821_F3,187_50.1,155_50.1,149_50.1,137_50.1,130_50.1,24_50.1,18_50.1

		String Sequence = "";
		while(ER.lookAhead() != '>' && ER.more())
			Sequence += ER.readLine();

		FastaSequence FS = new FastaSequence(InfoLine,Sequence,nrOfExperiments);
		this.add(FS);
	}


	public static void fixFastaFile(String inFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));
			ExtendedWriter EW = new ExtendedWriter(new FileWriter (inFile+".fixed"));

			int count = 0;
			while(ER.more()){

				if((char)ER.lookAhead() == '>'){
					String seqName = ER.readLine();
					String fixedSeqName = IOTools.fixFastaNames(seqName)+"_"+count;
					EW.println(fixedSeqName);
					count++;
				}
				else{
					EW.println(ER.readLine());
				}
			}
			EW.flush();
			EW.close();
			ER.close();
		}catch(Exception E){E.printStackTrace();}

	}

	public static void fixOneLine(String inFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));
			ExtendedWriter EW = new ExtendedWriter(new FileWriter (inFile+".fixed"));
			ExtendedWriter EW2 = new ExtendedWriter(new FileWriter (inFile+".size"));


			int count = 0;
			boolean first  = true;
			int length = 0;
			while(ER.more()){

				if((char)ER.lookAhead() == '>'){
					String seqName = ER.readLine();
					String fixedSeqName = IOTools.fixFastaNames(seqName)+"_"+count;
					if(count != 0 )EW.println();
					EW.println(fixedSeqName);
					count++;
					if(!first){
						EW2.println("\t"+length);
						length = 0;
					}
					else
						first = false;
					EW2.print(seqName);
				}
				else{
					String seq = ER.readLine();
					length += seq.length();
					EW.print(seq);
				}
			}
			EW2.println("\t"+length);
			EW.flush();
			EW.close();
			ER.close();
			EW2.flush();
			EW2.close();
		}catch(Exception E){E.printStackTrace();}

	}

	public static void fixArrowProblem(String inFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));
			ExtendedWriter EW = new ExtendedWriter(new FileWriter (inFile+".fixed"));

			int count = 0;
			boolean first  = true;
			int length = 0;
			while(ER.more()){

				if((char)ER.lookAhead() == '>'){
					if(count != 0 )EW.println();

					String seqName = ER.readLine();
					char[] Narray = seqName.toCharArray();
					int pointer = 0;
					while(Narray[pointer] == '>')pointer++;
					seqName = seqName.substring(pointer);
					EW.println(">"+seqName);
					count++;
					if(!first){
						length = 0;
					}
					else
						first = false;
				}
				else{
					String seq = ER.readLine();
					length += seq.length();
					EW.print(seq);
				}
			}
			EW.println();
			EW.flush();
			EW.close();
			ER.close();
		}catch(Exception E){E.printStackTrace();}

	}



	public static void extractSeqAboveLength(String inFile, String outFile, String WD, int length){
		try{
			FastaSequences.oneLine(inFile,WD+"/"+inFile+"_temp");
			ExtendedReader ER = new ExtendedReader(new FileReader(WD+"/"+inFile+"_temp"));
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile));
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					String name = ER.readLine();
					String seq = ER.readLine();
					if(seq.length()> length){
						EW.println(name);
						EW.println(seq);
					}
				}
			}
			IOTools.deleteFile(WD, inFile+"_temp");
			ER.close();
			EW.flush();
			EW.close();
		}catch(Exception E){E.printStackTrace();}
	}


	public static void oneLine(String inFile){
		oneLine(inFile,inFile+".oneLine" );
	}

	public static void oneLineFirst(String inFile){
		oneLineFirst(inFile,inFile+".oneLine" );
	}


	public static void oneLine(String inFile, String outFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile));


			int count = 0;
			boolean first  = true;
			int length = 0;
			while(ER.more()){

				if((char)ER.lookAhead() == '>'){
					String seqName = ER.readLine();
					if(count != 0 )EW.println();
					EW.println(seqName);
					count++;
					if(!first){
						length = 0;
					}
					else
						first = false;
				}
				else{
					String seq = ER.readLine();
					length += seq.length();
					EW.print(seq);
				}
			}
			EW.flush();
			EW.close();
			ER.close();
		}catch(Exception E){E.printStackTrace();}

	}

	public static void oneLineFirst(String inFile, String outFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outFile));


			int count = 0;
			boolean first  = true;
			int length = 0;
			while(ER.more()){

				if((char)ER.lookAhead() == '>'){
					String seqName = ER.readWord();
					ER.skipLine();
					if(count != 0 )EW.println();
					if(seqName.indexOf("/")== -1)EW.println(seqName);
					else {
						EW.println(seqName.substring(0,seqName.indexOf("/")));
					}
					count++;
					if(!first){
						length = 0;
					}
					else{
						first = false;

					}
				}
				else{
					String seq = ER.readLine();
					length += seq.length();
					EW.print(seq);
				}
			}
			EW.flush();
			EW.close();
			ER.close();
		}catch(Exception E){E.printStackTrace();}

	}



	public static void split(String inFile, int nrOfSequences){

		String[] splitFileNames = null;
		IOTools.mkDir("split");
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));


			int temp = nrOfSequences;
			int count = 0;
			int number = 0;
			while(ER.more()){
				ExtendedWriter EW = new ExtendedWriter(new FileWriter ("split/"+inFile+"_"+number+".fa"));
				while(ER.more() && count < nrOfSequences){

					if((char)ER.lookAhead() == '>'){
						String seqName = ER.readLine();
						EW.println(seqName);
						count++;
						while(ER.more() && (char)ER.lookAhead() != '>'){
							String seq = ER.readLine();
							EW.println(seq);
						}
					}
				}

				EW.flush();
				EW.close();
				nrOfSequences+= temp;
				number++;
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}

	}


	public static ArrayList <String> split(String inDir, String inFile, int nrOfSequences, String suffix){
		String inFileBase = IOTools.getFileBase(inFile, suffix);

		ArrayList <String> splitFileNames = new ArrayList<String>();
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inDir+"/"+inFile));

			int temp = nrOfSequences;
			int count = 0;
			int number = 0;
			while(ER.more()){

				ExtendedWriter EW = new ExtendedWriter(new FileWriter (inDir+"/"+inFileBase+"_"+number+"."+suffix));
				splitFileNames.add(inFileBase+"_"+number+"."+suffix);
				while(ER.more() && count < nrOfSequences){

					if((char)ER.lookAhead() == '>'){
						String seqName = ER.readLine();
						EW.println(seqName);
						count++;
						while(ER.more() && (char)ER.lookAhead() != '>'){
							String seq = ER.readLine();
							EW.println(seq);
						}
					}
				}

				EW.flush();
				EW.close();
				nrOfSequences+= temp;
				number++;
			}
			ER.close();
			return splitFileNames;
		}catch(Exception E){E.printStackTrace();}
		return null;
	}

	public static ArrayList <String> splitSize(String inDir, String inFile, int Size, String suffix){
		String inFileBase = IOTools.getFileBase(inFile, suffix);

		ArrayList <String> splitFileNames = new ArrayList<String>();
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inDir+"/"+inFile));

			int temp = Size;
			int length = 0;
			int number = 0;
			while(ER.more()){
				if(!IOTools.isDir(inDir+"/tmp"))
					IOTools.mkDir(inDir+"/tmp");
				ExtendedWriter EW = new ExtendedWriter(new FileWriter (inDir+"/tmp/"+inFileBase+"_"+number+"."+suffix));
				splitFileNames.add(inFileBase+"_"+number+"."+suffix);
				length = 0;
				while(ER.more() && length < Size){
					if((char)ER.lookAhead() == '>'){
						String seqName = ER.readLine();
						EW.println(seqName);
						while(ER.more() && (char)ER.lookAhead() != '>'){
							String seq = ER.readLine();
							EW.println(seq);
							length = length+seq.length();
						}
					}
				}
				EW.flush();
				EW.close();
				number++;
			}
			ER.close();
			return splitFileNames;
		}catch(Exception E){E.printStackTrace();}
		return null;
	}

	public static ArrayList <String> splitSizefileNames(String inDir, String inFile, int Size, String suffix){
		String inFileBase = IOTools.getFileBase(inFile, suffix);

		ArrayList <String> splitFileNames = new ArrayList<String>();
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inDir+"/"+inFile));

			int temp = Size;
			int length = 0;
			int number = 0;
			while(ER.more()){
				if(!IOTools.isDir(inDir+"/tmp"))
					IOTools.mkDir(inDir+"/tmp");
				splitFileNames.add(inFileBase+"_"+number+"."+suffix);
				length = 0;
				while(ER.more() && length < Size){
					if((char)ER.lookAhead() == '>'){
						String seqName = ER.readLine();
						while(ER.more() && (char)ER.lookAhead() != '>'){
							String seq = ER.readLine();
							length = length+seq.length();
						}
					}
				}
				number++;
			}
			ER.close();
			return splitFileNames;
		}catch(Exception E){E.printStackTrace();}
		return null;
	}

	public static void getLength(String inFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));
			ExtendedWriter EW2 = new ExtendedWriter(new FileWriter (inFile+".size"));


			int count = 0;
			boolean first  = true;
			int length = 0;
			while(ER.more()){

				if((char)ER.lookAhead() == '>'){
					String seqName = ER.readLine();
					seqName=seqName.split(" ")[0].substring(1);
					count++;
					if(!first){
						EW2.println("\t"+length);
						length = 0;
					}
					else
						first = false;
					EW2.print(seqName);
				}
				else{
					String seq = ER.readLine();
					length += seq.length();
				}
			}
			EW2.println("\t"+length);
			ER.close();
			EW2.flush();
			EW2.close();
		}catch(Exception E){E.printStackTrace();}

	}

	public static void getInfo(String inFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));
			ExtendedWriter EW2 = new ExtendedWriter(new FileWriter (inFile+".info"));


			int count = 0;
			boolean first  = true;
			int length = 0;
			int gc=0;
			char[] needles = "GC".toCharArray();
			while(ER.more()){

				if((char)ER.lookAhead() == '>'){
					String seqName = ER.readLine();
					seqName=seqName.split(" ")[0].substring(1);
					count++;
					if(!first){
						EW2.println("\t"+length+"\t"+((double)gc/(double)length));
						length = 0;
						gc=0;
					}
					else
						first = false;
					EW2.print(seqName);
				}
				else{
					String seq = ER.readLine().toUpperCase();
					gc += Functions.countOccurrences(seq, needles);
					length += seq.length();
				}
			}
			EW2.println("\t"+length+"\t"+((double)gc/(double)length));
			ER.close();
			EW2.flush();
			EW2.close();
		}catch(Exception E){E.printStackTrace();}

	}


	public static void RNA2DNA(String inFile){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));
			ExtendedWriter EW = new ExtendedWriter(new FileWriter (inFile+".DNA"));
			boolean first  = true;

			while(ER.more()){

				if((char)ER.lookAhead() == '>'){
					String seqName = ER.readLine();
					if(first)
						first = false;
					else
						EW.println();
					EW.println(seqName);
				}
				else{
					EW.print(RNAfunctions.RNA2DNA(ER.readLine()));
				}
			}
			EW.flush();
			EW.close();
			ER.close();
		}catch(Exception E){E.printStackTrace();}

	}


	public static void extractEnds(String inFile, int length){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));
			ExtendedWriter EW = new ExtendedWriter(new FileWriter (inFile+".ends"));
			boolean first  = true;
			ArrayList <String>sequence = new ArrayList<String>();
			String seqName = "";
			while(ER.more()){

				if((char)ER.lookAhead() == '>'){
					if(first)
						first = false;
					else{
						int i = 0;
						String Fend = "";
						while(i < sequence.size() && Fend.length() < length){
							Fend += sequence.get(i);
							i++;
						}
						if(Fend.length()> length)
							Fend = Fend.substring(0,length);
						EW.println(seqName+"_5end");
						EW.println(Fend);
						i = sequence.size()-1;
						String Tend = "";
						while(i > -1 && Tend.length() < length){
							Tend = sequence.get(i)+Tend;
							i--;
						}
						if(Tend.length()> length)
							Tend = Tend.substring(Tend.length()-length);
						EW.println(seqName+"_3end");
						EW.println(Tend);
					}
					seqName = ER.readLine();
					sequence = new ArrayList<String>();
				}
				else{
					sequence.add(ER.readLine());
				}
			}
			EW.flush();
			EW.close();
			ER.close();
		}catch(Exception E){E.printStackTrace();}

	}

	public static void padEnds(String inFile, String pad){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));
			ExtendedWriter EW = new ExtendedWriter(new FileWriter (inFile+".padded"));
			boolean first  = true;
			ArrayList <String>sequence = new ArrayList<String>();
			String seqName = "";
			//			while(ER.more()){
			//
			//				if((char)ER.lookAhead() == '>'){
			//					if(first)
			//						first = false;
			//					else{
			//						int i = 0;
			//						String Fend = "";
			//						while(i < sequence.size() && Fend.length() < length){
			//							Fend += sequence.get(i);
			//							i++;
			//						}
			//						if(Fend.length()> length)
			//							Fend = Fend.substring(0,length);
			//						EW.println(seqName+"_5end");
			//						EW.println(Fend);
			//						i = sequence.size()-1;
			//						String Tend = "";
			//						while(i > -1 && Tend.length() < length){
			//							Tend = sequence.get(i)+Tend;
			//							i--;
			//						}
			//						if(Tend.length()> length)
			//							Tend = Tend.substring(Tend.length()-length);
			//						EW.println(seqName+"_3end");
			//						EW.println(Tend);
			//					}
			//					seqName = ER.readLine();
			//					sequence = new ArrayList<String>();
			//				}
			//				else{
			//					sequence.add(ER.readLine());
			//				}
			//			}
			//			EW.flush();
			//			EW.close();
			//			ER.close();
		}catch(Exception E){E.printStackTrace();}

	}

	private void addFasta(ExtendedReader ER){
		String InfoLine = ER.readLine();


		int[] tempSeq = new int[10000000];
		int pointer = 0;

		while(ER.lookAhead() != '>' && ER.more()){
			String A = ER.readLine().trim();
			int[] subSequence = RNAfunctions.RNAString2Int(A);
			for(int i = 0; i < subSequence.length;i++){
				tempSeq[pointer] = subSequence[i];
				pointer++;
			}
			ER.more();

		}
		char[] A = new char[pointer];
		for(int i = 0; i < pointer; i++){
			A[i] = RNAfunctions.RNAInt2char(tempSeq[i]);
		}
		FastaSequence FS = new FastaSequence(InfoLine,new String(A));
		this.add(FS);
	}

	private static void getORFs(ExtendedReader ER,ExtendedWriter EW,ExtendedWriter EW2,ExtendedWriter EW3,ExtendedWriter EW4,ExtendedWriter EW5,GeneticCode GC){
		String InfoLine = ER.readLine();
		//>853_36_1821_F3,187_50.1,155_50.1,149_50.1,137_50.1,130_50.1,24_50.1,18_50.1

		int[] tempSeq = new int[10000000];
		int pointer = 0;

		while(ER.lookAhead() != '>' && ER.more()){
			String A = ER.readLine().trim();
			int[] subSequence = RNAfunctions.RNAString2Int(A);
			for(int i = 0; i < subSequence.length;i++){
				tempSeq[pointer] = subSequence[i];
				pointer++;
			}
			ER.more();

		}
		char[] A = new char[pointer];
		for(int i = 0; i < pointer; i++){
			A[i] = RNAfunctions.RNAInt2char(tempSeq[i]);
		}
		FastaSequence FS = null;
		if(InfoLine.indexOf(" ")==-1)
			FS = new FastaSequence(InfoLine,new String(A));
		else
			FS = new FastaSequence(InfoLine.substring(0,InfoLine.indexOf(" ")),new String(A));

		FS.findLongestORF(0,EW,EW2,EW3,EW4,EW5,GC);
	}



	private void addFastaName(ExtendedReader ER){
		String InfoLine = ER.readLine();
		int length = 0;
		while(ER.lookAhead() != '>' && ER.more()){
			String seq = ER.readLine();
			length+= seq.length();
			ER.more();
		}
		FastaSequence FS = new FastaSequence(InfoLine.substring(1));
		FS.length = length;
		this.add(FS);
	}

	public void printFasta(String dir, String outfile){
		try{
			ExtendedWriter EW= new ExtendedWriter(new FileWriter(dir+"/"+outfile));
			for(int i = 0; i < this.size();i++){
				this.get(i).printFasta(EW);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){
			E.printStackTrace();
		}
	}

	public void printFasta(ExtendedWriter EW){
		for(int i = 0; i < this.size();i++){
			this.get(i).printFasta(EW);
		}
	}


	public void printFastaRNA(String dir, String outfile){
		try{
			ExtendedWriter EW= new ExtendedWriter(new FileWriter(dir+"/"+outfile));
			for(int i = 0; i < this.size();i++){
				this.get(i).printFastaRNA(EW);
			}
			EW.flush();
			EW.close();
		}catch(Exception E){
			E.printStackTrace();
		}
	}

	public void printFastaSubset(String dir, String infile, String outfile, Hashtable<String,String> HT){
		try{
			ExtendedWriter EW= new ExtendedWriter(new FileWriter(dir+"/"+outfile));
			ExtendedReader ER= new ExtendedReader(new FileReader(dir+"/"+infile));

			this.printFastaSubset(ER, HT, EW);
			ER.close();
			EW.flush();
			EW.close();
		}catch(Exception E){
			E.printStackTrace();
		}
	}
	public void printFastaSubsetmiRNA(String dir, String infile, String outfile, Hashtable<String,String> HT){
		try{
			ExtendedWriter EW= new ExtendedWriter(new FileWriter(dir+"/"+outfile));
			ExtendedReader ER= new ExtendedReader(new FileReader(dir+"/"+infile));

			this.printFastaSubsetmiRNA(ER, HT, EW);
			ER.close();
			EW.flush();
			EW.close();
		}catch(Exception E){
			E.printStackTrace();
		}
	}





	private void printFastaSubset(ExtendedReader ER,Hashtable<String,String> HT, ExtendedWriter EW){
		while(ER.more()){
			String InfoLine = ER.readLine().substring(1);
			//System.out.println(InfoLine);
			String[] info = InfoLine.split("\\ ");
			info = info[0].split("/");
			if(!HT.containsKey(info[0])){
				while(ER.lookAhead() != '>' && ER.more()){
					ER.skipLine();
					ER.more();
				}
			}else{
				EW.println(">"+InfoLine);
				while(ER.lookAhead() != '>' && ER.more())
					EW.println( ER.readLine());
			}

		}

	}

	private void printFastaSubsetmiRNA(ExtendedReader ER,Hashtable<String,String> HT, ExtendedWriter EW){
		while(ER.more()){
			String InfoLine = ER.readLine().substring(1);
			System.out.println(InfoLine);
			String[] info = InfoLine.split("\\ ");
			String[] info2 = info[0].split("-");
			String name = "wrong";
			if(info2.length > 2)
				name = info2[0]+"-"+info2[1]+"-"+info2[2];

			if(!HT.containsKey(name)){
				while(ER.lookAhead() != '>' && ER.more()){
					ER.skipLine();
					ER.more();
				}
			}else{
				EW.println(">"+InfoLine);
				while(ER.lookAhead() != '>' && ER.more())
					EW.println( ER.readLine());
			}

		}

	}



	public static FastaSequence getFasta(ExtendedReader ER){
		String InfoLine = ER.readLine();
		String Sequence = "";
		while(ER.lookAhead() != '>' && ER.more())
			Sequence += ER.readLine();

		return new FastaSequence(InfoLine,Sequence);
	}



	public static void run(Hashtable<String,String> T){

		if(T.containsKey("-extra")){
			String extra = Functions.getValue(T, "-extra", "trinity");
			String dir =  Functions.getValue(T, "-d", ".");
			String inFile =  Functions.getValue(T, "-i", "fileName");
			String outFile =  Functions.getValue(T, "-o", inFile+".extra.fa");
			if(inFile.compareTo("fileName") != 0){
				FastaSequences FS = new FastaSequences(dir, inFile);
				FS.printSequences(dir, outFile,extra);
			}
		}
		if(T.containsKey("-subset")){
			String subset = Functions.getValue(T, "-subset", "subset.txt");
			String dir =  Functions.getValue(T, "-d", IOTools.getCurrentPath());
			String inFile =  Functions.getValue(T, "-i", "fileName");
			String outFile =  Functions.getValue(T, "-o", inFile+".subset.fa");
			System.out.println("writing sequences found in "+ inFile +" and "+subset+ "in file "+outFile);
			if(inFile.compareTo("fileName") != 0 ){
				Hashtable<String,String> HT = new Hashtable<String,String>(1000000);
				HT = Functions.getKeys(HT, subset);
				FastaSequences FS = new FastaSequences();
				if(T.containsKey("-miRNA"))
					FS.printFastaSubsetmiRNA(dir, inFile,outFile, HT);
				else
					FS.printFastaSubset(dir, inFile,outFile, HT);

			}
		}
		if(T.containsKey("-findMERs")){
			String mers = Functions.getValue(T, "-mers", "subset.txt");
			String dir =  Functions.getValue(T, "-d", IOTools.getCurrentPath());
			String inFile =  Functions.getValue(T, "-i", "searchSpace");
			String outFile =  Functions.getValue(T, "-o", inFile+"_"+mers+".matches");
			ExtendedWriter EW = ExtendedWriter.getFileWriter(dir+"/"+outFile);	
			System.out.println("finding mers from file "+mers+" in file "+inFile+" and printing them to "+outFile);
			if(inFile.compareTo("fileName") != 0 ){
				FastaSequences MERS = new FastaSequences(dir,mers);
				FastaSequences INFILE = new FastaSequences(dir,inFile);
				for(int i = 0; i < INFILE.getNames().length; i++){
					for(int j = 0; j < MERS.size(); j++){
						FastaSequence A = INFILE.get(i);
						FastaSequence B = MERS.get(j);
						ArrayList<Integer> locations = A.findMers(B);
						if(locations.size() > 0){
							for(int k = 0; k < locations.size();k++){
								if(locations.get(k) >-1)
									EW.println(A.getName()+"\t"+B.getName()+"\t+\t"+locations.get(k)+"\t"+(locations.get(k)+B.Sequence.length));
								else
									EW.println(A.getName()+"\t"+B.getName()+"\t-\t"+-locations.get(k)+"\t"+(-locations.get(k)-B.Sequence.length));
							}
						}
					}
				}
			}
			EW.flush();
			EW.close();
		}

	}


	public void getFastaFileNames(String dir, String fileName){
		getFastaFileNames(dir+"/"+fileName);
	}

	public void getFastaFileNames(String fileName){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			int count = 1000;
			while(ER.more()){

				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{	
					addFastaName(ER);
				}
				if(this.size() > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 1000;
				}
			}
			System.out.println("Finished total nr of sequences = " + this.size());
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}




	void readFastaFile(String dir, String fileName){
		readFastaFile(dir+"/"+fileName);
	}

	public void readFastaFile(String fileName){
		if(verbose)
			System.out.println("Verbose mode is on");
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			int count = 1000;
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{	
					addFasta(ER);
				}
				if(this.size() > count && verbose){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 1000;
				}
			}
			if(verbose)
				System.out.println("Finished total nr of sequences = " + this.size());
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}


	public static FastaSequences getFastaSequences(String inFile, String WD){
		FastaSequences FSs = new FastaSequences();
		FSs.parseAllFasta(inFile, WD);

		return FSs;
	}

	public void parseAllFasta(String inFile, String WD){
		try{
			FastaSequences.oneLineFirst(WD+"/"+inFile,WD+"/temp.fa");
			ExtendedReader ER = new ExtendedReader(new FileReader(WD+"/temp.fa"));
			int count = 1000;
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					addAllFasta(ER);
				}
				if(this.size() > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 1000;
				}
			}
			System.out.println("Finished total nr of sequences = " + this.size());
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}


	private void addAllFasta(ExtendedReader EW){
		while(EW.more()){
			String infoLine = EW.readLine();
			String[] info = infoLine.split("\\ ");

			String hitName = info[0].substring(1);
			String Sequence = EW.readLine();
			FastaSequence A = new FastaSequence(hitName, Sequence);
			this.add(A);
		}
	}

	public ArrayList<String> getAllSequenceNames(){
		ArrayList <String> ContigNames = new ArrayList<String>();

		for(int i = 0; i < this.size(); i++ ){
			ContigNames.add(this.get(i).Name);
		}
		return ContigNames;
	}


	void readFastaFile(String dir, String fileName, int nrOfExperiments){

		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (dir+"/"+fileName));
			int count = 100000;
			while(ER.more()){

				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{
					addFasta(ER, nrOfExperiments);
				}
				if(this.size() > count){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 100000;
				}
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}


	public  static void getORFs(String fileName){
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader(fileName));
			int split = fileName.lastIndexOf(".");
			String fileBase = fileName.substring(0,split);
			String suffix = fileName.substring(split+1);
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(fileBase+".info"));
			ExtendedWriter EW2 = new ExtendedWriter(new FileWriter(fileBase+".ORFs."+suffix));
			ExtendedWriter EW3 = new ExtendedWriter(new FileWriter(fileBase+".523prime."+suffix));
			ExtendedWriter EW4 = new ExtendedWriter(new FileWriter(fileBase+".peptide."+suffix));
			ExtendedWriter EW5 = new ExtendedWriter(new FileWriter(fileBase+".subset."+suffix));

			EW.println("Name\tDir\tORFlength\tSeqLength\tstart\tstop\tkind\tGC_count");

			GeneticCode GC = new GeneticCode();
			int count = 1000;
			int pointer = 0;
			while(ER.more()){
				if((char)ER.lookAhead() == '#'){
					ER.skipLine();
				}
				else{	
					getORFs(ER,EW,EW2,EW3,EW4,EW5,GC);
					pointer++;
				}
				if(pointer > count ){
					System.out.println("nr Of Sequences Read" + count);
					count= count + 1000;
				}
			}
			System.out.println("Finished total nr of sequences = " + pointer);
			EW.flush();
			EW.close();
			EW2.flush();
			EW2.close();
			EW3.flush();
			EW3.close();
			EW4.flush();
			EW4.close();
			ER.close();
		}catch(Exception E){E.printStackTrace();}
	}




	private void printInfo(String cfastaDir , String cfastaFile, String compareFastaFile, String compareDir, int nrOfSequences, int nrOfcomparingSequences, int nrOfuniqueSequences, String outDir,String outFile){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile+  "  .......");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile));
			EW.println("# original file =  "+ cfastaDir +"/"+cfastaFile );
			EW.println("# compare file = "+ compareDir +"/"+compareFastaFile );
			EW.println("# Nr of sequences in original file :"+nrOfSequences );
			EW.println("# Nr of sequences in compare file :"+nrOfcomparingSequences );
			EW.println("# Nr of unique sequences in original file :"+nrOfuniqueSequences );
			EW.println("# Nr of unique sequences in compare file :"+(nrOfcomparingSequences - (nrOfSequences - nrOfuniqueSequences)));
			EW.println("# Nr of same seqences in both files :"+(nrOfSequences - nrOfuniqueSequences) );

			for(int i = 0; i < this.size(); i++){
				this.get(i).printFasta(EW);
			}
			EW.flush();
			EW.close();

		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}


	private void printNrOfHits(int nrOfCfastaFiles,String outDir, String outFile){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile+  "  .......");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile));
			EW.println("# Nr of cFastaSequences:"+nrOfCfastaFiles);

			for(int i = 0; i < this.size(); i++){
				this.get(i).printNrOfHits(EW);
			}
			EW.flush();
			EW.close();

		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}

	private void printNrOfHits(String outDir, String outFile, String[] experiments){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile+  "  .......");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile));
			EW.print("Name\tlength");
			for(int i = 0; i < experiments.length; i++){
				EW.print("\t"+experiments[i]);
			}
			EW.println();
			for(int i = 0; i < this.size(); i++){
				this.get(i).printNrOfHits(EW);
			}
			EW.flush();
			EW.close();

		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}


	private void printSequences(ExtendedWriter EW){
		for(int i = 0; i < this.size(); i++){
			this.get(i).printFasta(EW);
		}
	}

	public void padSequences(String Sequence){
		for(int i = 0; i < this.size(); i++){
			this.get(i).padSequenceBothEnds(Sequence);
		}
	}

	private void printSequences(ExtendedWriter EW, String Extra){
		for(int i = 0; i < this.size(); i++){
			this.get(i).printFasta(EW, Extra+"_"+i);
		}
	}

	private void printSequences(String outDir,String outFile){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile+"  .......");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile));
			printSequences(EW);
			EW.flush();
			EW.close();

		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}

	private void printSequences(String outDir,String outFile,String Extra){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile+"  .......");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile));
			printSequences(EW, Extra);
			EW.flush();
			EW.close();

		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}


	private void printSingle(String outDir,String outFile){
		System.out.print("Printing sequences to file "+outDir+"/"+outFile+  ".single  .......");
		try{
			ExtendedWriter EW = new ExtendedWriter(new FileWriter(outDir+"/"+outFile+".single"));
			EW.println("# fasta sequence "+outDir +"/"+ outFile +" with one unique hit :");
			EW.println("#Total number of sequences with one hit: "+  this.size());
			printSequences(EW);
			EW.flush();
			EW.close();

		}catch(Exception E){E.printStackTrace();}
		System.out.println("finished");
	}


	public void printSequence(String Name, int start, int stop){

		for(int i = 0; i < this.size(); i++){
			if(this.get(i).Name.indexOf(Name) > 0){
				System.out.println(this.get(i).Name+"_"+start+"->"+stop);
				System.out.println(this.get(i).getSequence(start, stop));
			}
		}
	}

	public void printSequenceSurr(String Name, int start, int stop, int surr,String extra){

		for(int i = 0; i < this.size(); i++){
			if(this.get(i).Name.indexOf(Name) > 0){
				System.out.println(this.get(i).Name+" ("+start+"->"+stop+") "+extra);
				System.out.println(this.get(i).getSequenceSurr(start, stop,surr));
			}
		}
	}


	/*	private final int mRNA = 1;
	private final int ncRNA = 2;
	private final int intergenic = 3;
	private final int antisense = 4;
	 */

	public static void printDistribution(int[] Distribution, int start, int step){
		for(int i = 0; i< Distribution.length; i++)
			System.out.println((start+i*step)+"\t"+Distribution[i]);
	}



	public void setName(String name) {
		Name = name;
	}

	public String getName() {
		return Name;
	}


	public int[] findGene(String GeneName){
		int[] positions = null;
		if(this.size() == 0) return null;
		for(int i = 0; i < this.size();i++){
			if(this.get(i).getName().indexOf(GeneName) > -1)
				positions = Functions.addInt(positions,i);
		}
		if(positions == null){
			for(int i = 0; i < this.size();i++){
				if(GeneName.indexOf(this.get(i).getName().split("\t")[0]) != -1){
					positions = Functions.addInt(positions,i);
					return positions;
				}
			}
		}
		return positions;
	}

	public int[] findGene(int[] Sequence){
		int[] positions = null;
		if(this.size() == 0) return null;
		for(int i = 0; i < this.size();i++){
			if(RNAfunctions.compareSequences(Sequence, this.get(i).Sequence))
				positions = Functions.addInt(positions,i);
		}
		return positions;
	}



	public static FastaSequences getUniqueNames(FastaSequences FS1, FastaSequences FS2){
		FastaSequences Unique = new FastaSequences();
		for(int i = 0; i <  FS1.size();i++){
			int[] location = FS2.findGene(FS1.get(i).Name);
			if(location == null)
				Unique.add(FS1.get(i)); 

		}
		return Unique;
	}


	public static FastaSequences getUniqueSequences(FastaSequences FS1, FastaSequences FS2){
		FastaSequences Unique = new FastaSequences();
		for(int i = 0; i <  FS1.size();i++){
			int[] location = FS2.findGene(FS1.get(i).Sequence);
			if(location == null)
				Unique.add(FS1.get(i));
		}
		return Unique;
	}

	public static FastaSequences merge(FastaSequences FS1, FastaSequences FS2){
		FastaSequences both = new FastaSequences();
		for(int i = 0; i < FS2.size(); i++)
			both.add(FS2.get(i));
		for(int i = 0; i <  FS1.size();i++){
			both.add(FS1.get(i));
		}
		return both;
	}


	public void writeFastaFile(String Dir,String FileName){
		try{
			ExtendedWriter ew = new ExtendedWriter(new FileWriter(Dir+"/"+FileName));
			for(int i = 0 ; i < this.size(); i++){
				ew.println(">"+ this.get(i).getName());
				ew.println(RNAfunctions.RNAInt2String(this.get(i).getSequence()));
			}
			ew.flush();
			ew.close();
		}
		catch(Exception E){
			E.printStackTrace();
		}
	}

	public void writeFastaFileDNA(String Dir,String FileName){
		try{
			ExtendedWriter ew = new ExtendedWriter(new FileWriter(Dir+"/"+FileName));
			for(int i = 0 ; i < this.size(); i++){
				ew.println(">"+ this.get(i).getName());
				ew.println(RNAfunctions.DNAInt2String(this.get(i).getSequence()));
			}
			ew.flush();
			ew.close();
		}
		catch(Exception E){
			E.printStackTrace();
		}
	}

	public int[] getLengths(){
		int[] lengths = new int[this.size()];
		for(int i = 0; i < this.size(); i++){
			lengths[i] = this.get(i).getSequence().length;
		}
		return lengths;
	}

	public String[] getNames(){
		String[] Names = new String[this.size()];
		for(int i = 0; i < this.size(); i++){
			Names[i] = this.get(i).Name;
		}
		return Names;
	}




	public static ArrayList<Integer> getLengths(String inFile){

		ArrayList<Integer> Lengths = new ArrayList<Integer>();
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));


			boolean first  = true;
			int length = 0;
			while(ER.more()){

				if((char)ER.lookAhead() == '>'){
					String seqName = ER.readLine();
					seqName=seqName.substring(1);
					if(!first){
						Lengths.add(length);
						length = 0;
					}
					else
						first = false;
				}
				else{
					String seq = ER.readLine();
					length += seq.length();
				}
			}
			Lengths.add(length);

			ER.close();
		}catch(Exception E){E.printStackTrace();}
		return Lengths;

	}

	public static ArrayList<String> getNames(String inFile){

		ArrayList<String> Names = new ArrayList<String>();
		try{
			ExtendedReader ER = new ExtendedReader(new FileReader (inFile));


			while(ER.more()){
				if((char)ER.lookAhead() == '>'){
					String seqName = ER.readLine();
					seqName=seqName.substring(1);
					Names.add(seqName);
				}
				else{
					ER.readLine();
				}
			}
			ER.close();
		}catch(Exception E){E.printStackTrace();}
		return Names;

	}



	public int[] print5end(ExtendedWriter EW){
		int[] lengths = new int[this.size()];
		for(int i = 0; i < this.size(); i++){
			int[] sequence = this.get(i).getSequence();
			int[] subsequence = Functions.getSubarray(sequence, 0, 5);
			String name = this.get(i).Name;
			if(name.indexOf("*") < 0){
				EW.println(name+"_5end");
				EW.println(RNAfunctions.RNAInt2String(subsequence));
			}
		}
		return lengths;
	}


	public int[] print3end(ExtendedWriter EW){
		int[] lengths = new int[this.size()];
		for(int i = 0; i < this.size(); i++){
			int[] sequence = this.get(i).getSequence();
			int[] subsequence = Functions.getSubarray(sequence, sequence.length-5);
			String name = this.get(i).Name;
			if(name.indexOf("*") < 0){
				EW.println(name+"_3end");
				EW.println(RNAfunctions.DNAInt2String(subsequence));
			}
		}
		return lengths;
	}


	public int[] getSequence(String Name,int start, int stop){
		int[] position = findGene(Name);
		if(position != null){
			int[] fastaSeq = this.get(position[0]).getSequence();
			if(start < stop)
				return Functions.getSubarray(fastaSeq, start-1, stop);
			else{
				int[] revComp = Functions.getSubarray(fastaSeq, stop-1, start);
				return RNAfunctions.getReverseComplement(revComp);
			}
		}
		System.out.println(Name);
		return null;
	}

	public int[] getSequence(String Name,int start, int stop, int upstream, int downstream){
		int[] position = findGene(Name);
		if(position != null);
		int[] fastaSeq = this.get(position[0]).getSequence();
		if(start < stop){
			if(start - upstream > 0)
				start = start - upstream;
			else
				start = 1;
			if(stop + downstream < fastaSeq.length )
				stop = stop + downstream;
			else
				stop = fastaSeq.length;
			return getSequence(Name,start,stop);
		}
		else{
			if(stop - downstream > 0)
				stop = stop - downstream;
			else
				stop = 1;
			if(start + upstream < fastaSeq.length )
				start = start + upstream;
			else
				start = fastaSeq.length;
			return getSequence(Name,start,stop);
		}
	}


}




