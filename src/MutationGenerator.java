
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.StringTokenizer;

//import org.biojava3.core.sequence.ProteinSequence;
//import org.biojava3.core.sequence.io.FastaReaderHelper;


public class MutationGenerator {

	public static void main(String[] args) throws Exception{

		AAcid A = new AAcid('A', 1.8, 0.0);		//Ala
		AAcid R = new AAcid('R', -4.5 ,52.0);	//Arg
		AAcid N = new AAcid('N', -3.5, 3.380);	//Asn
		AAcid D = new AAcid('D', -3.5, 49.7);	//Asp
		AAcid C = new AAcid('C', 2.5, 1.480);	//Cys
		AAcid Q = new AAcid('Q', -3.5, 3.530);	//Gln
		AAcid E = new AAcid('E', -3.5, 49.90);	//Glu
		AAcid G = new AAcid('G', -0.4, 0.0);	//Gly
		AAcid H = new AAcid('H', -3.2, 51.600);	//His
		AAcid I = new AAcid('I', 4.5, 0.130);	//Ile
		AAcid L = new AAcid('L', 3.8, 0.130);	//Leu
		AAcid K = new AAcid('K', -3.9, 49.500);	//Lys
		AAcid M = new AAcid('M', 1.9, 1.430);	//Met
		AAcid F = new AAcid('F', 2.8, 0.350);	//Phe
		AAcid P = new AAcid('P', -1.6, 1.580);	//Pro
		AAcid S = new AAcid('S', -0.8, 1.670);	//Ser
		AAcid T = new AAcid('T', -0.7, 1.660);	//Thr		
		AAcid W = new AAcid('W', -0.9, 2.100);	//Trp
		AAcid Y = new AAcid('Y', -1.3, 1.610);	//Tyr
		AAcid V = new AAcid('V', 4.2, 0.130);	//Val

		//source: http://www.clcsupport.com/clcgenomicsworkbench/current/index.php?manual=Hydrophobicity_scales.html

		if(args.length < 4)
			printError();

		
		Scanner sc = null;
		try {
			sc = new Scanner(new File("hydroph.csv"));
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		String line = sc.nextLine();
		int counter = setHydrophobicityTable(args);	

		ArrayList<Integer> hydrophobicRegions = new ArrayList<Integer>();
		ArrayList<Integer> hydrogenBonds =  new ArrayList<Integer>();
		ArrayList<AAcid> aacids;
		ArrayList<AAcid> chain= new ArrayList<AAcid>();;
		aacids = new ArrayList<AAcid>();
		aacids.add(A);aacids.add(R);aacids.add(N);aacids.add(D);aacids.add(C);aacids.add(Q);aacids.add(E);aacids.add(G);
		aacids.add(H);aacids.add(I);aacids.add(L);aacids.add(K);aacids.add(M);aacids.add(F);aacids.add(P);aacids.add(S);
		aacids.add(T);aacids.add(W);aacids.add(Y);aacids.add(V);

		double hydrophr[] = new double[aacids.size()];  //hydrophobicity replace

		hydrophr = setNewHydrophobicity(hydrophr, sc, counter, aacids);

		/*
		for(int i = 0; i < aacids.size(); i++){
			System.out.println("aminoacid " + aacids.get(i).name() + " = " + aacids.get(i).hydrophobicity());
		}
		 */

		String input = "";
		System.out.println(args.length + " arguments entered");
		
		/*for(int i = 0; i < args.length; i++)
			System.out.println("argument " + i + " is " + args[i]);*/
			
		for (int i = 0; i < 2; i++){
			input = input + args[i];
			if(i !=1)
				input = input + ";";
		}

		StringTokenizer st = new StringTokenizer(input, ";"); 
		if( st.countTokens() != 2){
			printError();
		}
		else{
			//System.out.println("your input was : " + input);

			String hydroph = st.nextToken();
			String hydrob = st.nextToken();
			//	System.out.println("arg 1 = " + hydroph + " arg2 = " + hydrob);

			hydrophobicRegions = tokenize(hydroph);
			System.out.println("hydrophobic in the areas = " + hydrophobicRegions.toString());
			hydrogenBonds = tokenize(hydrob);
			System.out.println("hydrogen bonds in the areas = " + hydrogenBonds.toString());

		}
		
		String fasta = "";
		System.out.println("reading from file: " + args[3]);
		for(int i = 0; i < args.length; i++){
			if(args[i].equalsIgnoreCase("-chain") && ((i+1) <args.length)){
				fasta = args[i+1];
				break;
			}
			else if(args[i].contains(".pdb")){
				fasta = pdbToString(args[i]);
				break;
			}
		}
		if (fasta.equalsIgnoreCase(""))
				fasta = pdbToString(args[3]);
		//String fasta = args[2];  //the new input is the aminoacid chain
		System.out.println("the chain is " + fasta);

		chain = stringToAAcid(fasta, chain, aacids);  //convert the string to chain of AAcid objects

		double tendency1[] = new double[fasta.length()]; //hydrophobicity
		double tendency2[] = new double[fasta.length()]; //polarity

		for(int i = 0; i < fasta.length(); i++){
			if(hydrophobicRegions.contains(i) )
				tendency1[i] = 12.3 - chain.get(i).hydrophobicity();
			if(hydrogenBonds.contains(i)){
				tendency2[i] = 52.0 - chain.get(i).polarity();
			}
		}

		//there should be a formula that finds the balance between hydrophobicity and hydrogen bonding

		double tendency[] = new double[fasta.length()];
		for (int i = 0; i < fasta.length(); i++ )
			tendency[i] = tendency1[i] + tendency2[i];

		//determine the greatest tendency
		int greatestIndex = 0;
		double greatest = 0.0;
		for (int i = 0; i < fasta.length(); i++ ){
			if(tendency[i] > greatest){
				greatest = tendency[i];
				greatestIndex = i; 
			}
		}

		System.out.println("The change should be at " + greatestIndex + " the aminoacid: " + chain.get(greatestIndex).name() );


	}

	public static void printError(){
		System.out.println("Please enter the desired shape in the form of: hydrophobic regions; hydrogen bonded regions Hydrophobic_table_to_be_used AMINOACID_CHAIN");
		System.out.println("For example: 26-31,34 23,25 KD 2BEG.pdb ");
		System.out.println("Another example: 26-31,34 23,25 KD -chain DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA ");
		System.exit(1);
	}

	public static int setHydrophobicityTable(String args[]){
		int counter = 0;
		int argAt = 2;
		if (args[argAt].equalsIgnoreCase("KD")){
			counter = 2;
		}
		else if (args[argAt].equalsIgnoreCase("HW")){
			counter = 3;
		}
		else if (args[argAt].equalsIgnoreCase("Cor")){
			counter = 4;
		}
		else if (args[argAt].equalsIgnoreCase("Eis")){
			counter = 5;
		}
		else if (args[argAt].equalsIgnoreCase("Ros")){
			counter = 6;
		}
		else if (args[argAt].equalsIgnoreCase("Jan")){
			counter = 7;
		}
		else if (args[argAt].equalsIgnoreCase("EGES")){
			counter = 8;
		}
		return counter;
	}

	public static ArrayList<Integer> tokenize(String toke){
		ArrayList<Integer> list = new ArrayList<Integer>();

		StringTokenizer st = new StringTokenizer(toke, ",");  
		int tokeLength = st.countTokens();
		String tokearray[] = new String[tokeLength];
		for(int i= 0; i < tokeLength; i++){
			tokearray[i] = st.nextToken();
			//System.out.println(" tokearray [" + i + "] = " + tokearray[i]);
			if(tokearray[i].toString().contains("-")){ //if there is a range
				StringTokenizer st2 = new StringTokenizer(tokearray[i], "-");
		
				if (st2.countTokens() !=2)
					printError();
				int a=0,b=0;
				try{
					a = Integer.parseInt(st2.nextToken());
					b = Integer.parseInt(st2.nextToken());
					//	System.out.println(" a = " + a + " b = " + b);
				}catch(Exception e){printError();}
				if (a>b){
					int c = a;
					a = b;
					b = c;
				}
				for(int j = a; j <= b; j++ ){
					list.add(j);
				}
			}
			else{			//if there is not a range, only number
				try{
					list.add(Integer.parseInt(tokearray[i]));
				}catch(Exception e){}
			}
		}
		//System.out.println("arraylist = " + list.toString());




		return list;
	}

	public static double[] setNewHydrophobicity(double[] hydrophr, Scanner sc, int counter, ArrayList<AAcid> aacids){
		String line;
		int a = 0;
		while (sc.hasNextLine()){
			line = sc.nextLine();
			StringTokenizer st = new StringTokenizer(line, ",");
			st.nextToken();st.nextToken();
			for (int i = 1; i < counter; i++){
				hydrophr[a] = Double.parseDouble(st.nextToken());
				aacids.get(a).hphob = hydrophr[a];
			}
			a++;
		}
		sc.close();
	//	for(int i = 0; i < hydrophr.length; i++)
		//	System.out.println("entry " + i + " is " + hydrophr[i]);
		
		return hydrophr;

	}

	public static ArrayList<AAcid> stringToAAcid(String fasta, ArrayList<AAcid> chain, ArrayList<AAcid> aacids){
		for(int i = 0; i < fasta.length(); i++){     //convert the fasta seq to AAcid object arraylist
			for(int j = 0; j < aacids.size(); j++){
				Character ch = fasta.charAt(i);

				if((ch - aacids.get(j).name()) == 0){
					//	System.out.println("the character at " + i + " is " + ch);
					chain.add(aacids.get(j));
					break;
				}
			}
		}
		return chain;
	}

	public static String pdbToString(String filename){
		String pdbString = "";
		
		//reading from atom doesnt work because it has missing sequence. have to read from seqres. atom part has information containing the structure.
		
		try {
			Scanner sc2 = new Scanner(new File(filename));
			String line;
			int currentLine = 0;
			while(sc2.hasNextLine()){
				line = sc2.nextLine();
				if (line.regionMatches(true, 0, "SEQRES", 0, 6)){
					StringTokenizer tk = new StringTokenizer(line, " ");
					tk.nextToken();
					int lineNo = Integer.parseInt(tk.nextToken());
					if (lineNo < currentLine)
						break;
					currentLine++;
					//System.out.println(lineNo);
					tk.nextToken();
					String aStr = tk.nextToken();
				//	System.out.println(aStr);
					while(tk.hasMoreTokens()){
						pdbString = pdbString + CodeToLetter(tk.nextToken());
						//System.out.println(pdbString);
					}
				}
			}
			
		}catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}
		
		
			return pdbString;

		
	/*	try {
			Scanner sc2 = new Scanner(new File(filename));
			String line;
			String ch[] = new String[100];
			for(int i = 0; i < ch.length; i++)
				ch[i] = "";

			int[][] aNoS = new int[100][100];   //to be changed. first array is for different chains. second for the length of each
			int currentNo = 0;
			int a = 0;
			int b = 0;
			while(sc2.hasNextLine()){
				line = sc2.nextLine();
				if (line.regionMatches(true, 0, "ATOM", 0, 4)){
					StringTokenizer tk = new StringTokenizer(line, " ");
					tk.nextToken();tk.nextToken();tk.nextToken();
					String aStr = tk.nextToken();tk.nextToken();
					int aNo = Integer.parseInt(tk.nextToken());
					if (aNo == currentNo){
						continue;
					}
					else if(aNo > currentNo){
						currentNo = aNo;
						b++;
					}
					else{
						//this is a new string
						System.out.println("chain is " + (pdbString = ch[a]) );
						
						//a++;
						//currentNo = 0;
						//b = 0;
						//for now ignore the rest
						
						break;
					}
					System.out.println("a is " + a + " b is " + b);
					ch[a] = ch[a] + CodeToLetter(aStr);
					aNoS[a][b] = aNo;
					System.out.println("at aminoacid #" + aNo + " there is " + aStr + " first letter " + CodeToLetter(aStr) + " line " + line);
				}
			}

			sc2.close();
			int aacidLocations[] = new int[ch[a].length()];
			int c = 0;
			for (int i = 0; i < aNoS[a].length; i++){ //this will only be used if there is an issue with AminoAcids Not matching with the right locations
				if (aNoS[a][i] != 0){
					//System.out.print((aacidLocations[c] = aNoS[a][i]) + "-");
					c++;
				}
			}

			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
*/

	}

	public static char CodeToLetter(String code){
		char letter = ' ';

		if(code.equalsIgnoreCase("ALA")){
			letter = 'A';
		}else if(code.equalsIgnoreCase("ARG")){
			letter = 'R';
		}else if(code.equalsIgnoreCase("ASN")){
			letter = 'N';
		}else if(code.equalsIgnoreCase("ASP")){
			letter = 'D';
		}else if(code.equalsIgnoreCase("CYS")){
			letter = 'C';
		}else if(code.equalsIgnoreCase("GLN")){
			letter = 'Q';
		}else if(code.equalsIgnoreCase("GLU")){
			letter = 'E';
		}else if(code.equalsIgnoreCase("GLY")){
			letter = 'G';
		}else if(code.equalsIgnoreCase("HIS")){
			letter = 'H';
		}else if(code.equalsIgnoreCase("ILE")){
			letter = 'I';
		}else if(code.equalsIgnoreCase("LEU")){
			letter = 'L';
		}else if(code.equalsIgnoreCase("LYS")){
			letter = 'K';
		}else if(code.equalsIgnoreCase("MET")){
			letter = 'M';
		}else if(code.equalsIgnoreCase("PHE")){
			letter = 'F';
		}else if(code.equalsIgnoreCase("PRO")){
			letter = 'P';
		}else if(code.equalsIgnoreCase("SER")){
			letter = 'S';
		}else if(code.equalsIgnoreCase("THR")){
			letter = 'T';
		}else if(code.equalsIgnoreCase("TRP")){
			letter = 'W';
		}else if(code.equalsIgnoreCase("TYR")){
			letter = 'Y';
		}else if(code.equalsIgnoreCase("VAL")){
			letter = 'V';
		}

		return letter;
	}


}


