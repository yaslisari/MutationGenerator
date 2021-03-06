
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.StringTokenizer;

//Mohamed testing github here// 
//end of testing//


public class MutationGenerator {

	static int aacidLocations[];
	static double scalar1;
	static double scalar2;
	static int possibilities = 4;

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

		scalar1 = 1.0;
		scalar2 = 1.0;

		if(args.length < 4)
			printError();

		int counter = setHydrophobicityTable(args);	

		ArrayList<AAcid> aacids;
		ArrayList<AAcid> chain= new ArrayList<AAcid>();;
		aacids = new ArrayList<AAcid>();
		aacids.add(A);aacids.add(R);aacids.add(N);aacids.add(D);aacids.add(C);aacids.add(Q);aacids.add(E);aacids.add(G);
		aacids.add(H);aacids.add(I);aacids.add(L);aacids.add(K);aacids.add(M);aacids.add(F);aacids.add(P);aacids.add(S);
		aacids.add(T);aacids.add(W);aacids.add(Y);aacids.add(V);
		ArrayList[] arguments;

		double hydrophr[] = new double[aacids.size()];  //hydrophobicity replace

		hydrophr = setNewHydrophobicity(hydrophr, counter, aacids);
		char[] sortedHydrophobics = sortHydrophobics(aacids);
		char[] sortedPolarity = sortPolarity(aacids);

		arguments = manageArgs(args);			//parse arguments
		ArrayList<Integer> hydrophobicRegions = arguments[0];
		ArrayList<Integer> hydrogenBonds = arguments[1];
		ArrayList<Integer> hydrophilicRegions = arguments[2];
		ArrayList<Integer> noHydrogenBonds = arguments[3];

		String fasta = getFastaChain(args);		//get the chain from file or input

		chain = stringToAAcid(fasta, chain, aacids);  //convert the string to chain of AAcid objects

		//sort for the tendencies
		ArrayList<ChainAAcid> sortedTendencies = new ArrayList<ChainAAcid>();

		sortedTendencies = sort(sortedTendencies, fasta, chain, hydrophobicRegions, hydrogenBonds, hydrophilicRegions, noHydrogenBonds, aacids);


		char predicted[] = predictChange(sortedTendencies, sortedHydrophobics, sortedPolarity);
		//convert tendencies to string in the format of [(4,I) (36,S)]
		String output = tendencyToString(sortedTendencies, predicted);

	//	String path ="c:\\Windows\\System32\\cmd.exe";
	//	String input1 = "dir\n exit\n";

	//	System.out.println(sendToCommandline(path, input1));


	}

	private static ArrayList[] manageArgs(String[] args) {
		String input = "";
		System.out.println(args.length + " arguments entered");
		ArrayList<Integer> hydrophobicRegions = new ArrayList<Integer>();
		ArrayList<Integer> hydrogenBonds =  new ArrayList<Integer>();
		ArrayList<Integer> hydrophilicRegions = new ArrayList<Integer>();
		ArrayList<Integer> noHydrogenBonds = new ArrayList<Integer>();
		ArrayList arguments[] = new ArrayList[4];

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
			System.out.println("arg 1 = " + hydroph + " arg2 = " + hydrob);

			hydrophobicRegions = tokenize(hydroph);

			//take out hydrophillic regions
			for	(int i = 0; i < hydrophobicRegions.size(); i++){
				if(hydrophobicRegions.get(i)<0){
					hydrophilicRegions.add(hydrophobicRegions.get(i)*(-1));
					hydrophobicRegions.remove(i);
				}
			}

			System.out.println("hydrophobic in the areas = " + hydrophobicRegions.toString());
			System.out.println("hydrophilic in the areas = " + hydrophilicRegions.toString());

			hydrogenBonds = tokenize(hydrob);

			//take out no-hydrogen bond regions
			for	(int i = 0; i < hydrogenBonds.size(); i++){
				if(hydrogenBonds.get(i)<0){
					noHydrogenBonds.add(hydrogenBonds.get(i)*(-1));
					hydrogenBonds.remove(i);
				}
			}

			System.out.println("hydrogen bonds in the areas = " + hydrogenBonds.toString());
			System.out.println("no hydrogen bonds in the areas = " + noHydrogenBonds.toString());
			arguments[0] = hydrophobicRegions;
			arguments[1] = hydrogenBonds;
			arguments[2] = hydrophilicRegions;
			arguments[3] = noHydrogenBonds;
		}
		return arguments;
	}

	private static String getFastaChain(String[] args) {
		String fasta = "";
		//String fasta = args[2];  //the new input is the aminoacid chain

		//read chain from input or file
		for(int i = 0; i < args.length; i++){
			if(args[i].equalsIgnoreCase("-chain") && ((i+1) <args.length)){//if chain is entered manually
				fasta = args[i+1];
				break;
			}
			else if(args[i].contains(".pdb")){ //if chain is in the PDB file
				System.out.println("reading from file: " + args[3]);
				fasta = pdbToString(args[i]);
				break;
			}
		}

		System.out.println("the chain to be used is " + fasta);

		//	if (fasta.equalsIgnoreCase(""))  // burasinin islevi ne cozemedim.. eger gerekirse geri al
		//		fasta = pdbToString(args[3]);
		return fasta;
	}

	public static void printError(){
		System.out.println("Please enter the desired shape in the form of: hydrophobic regions; hydrogen bonded regions Hydrophobic_table_to_be_used AMINOACID_CHAIN");
		System.out.println("For example: 26:31,34,-35 23,25,-24 KD 2BEG.pdb ");
		System.out.println("Another example: 26:31,34,-35 23,25,-24 KD -chain DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA ");
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
			if(tokearray[i].toString().contains(":")){ //if there is a range
				StringTokenizer st2 = new StringTokenizer(tokearray[i], ":");

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

	public static double[] setNewHydrophobicity(double[] hydrophr, int counter, ArrayList<AAcid> aacids){
		Scanner sc = null;
		//read hydrophobicity tables
		try {
			sc = new Scanner(new File("hydroph.csv"));
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		String line = sc.nextLine();
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


		/*
	for(int i = 0; i < aacids.size(); i++){
		System.out.println("aminoacid " + aacids.get(i).name() + " = " + aacids.get(i).hydrophobicity());
	}
		 */

		return hydrophr;

	}

	public static ArrayList<AAcid> stringToAAcid(String fasta, ArrayList<AAcid> chain, ArrayList<AAcid> aacids){
		for(int i = 0; i < fasta.length(); i++){     //convert the fasta seq to AAcid object arraylist
			for(int j = 0; j < aacids.size(); j++){
				Character ch = fasta.charAt(i);

				if((ch - aacids.get(j).name()) == 0){
					//	System.out.println("the character at " + i + " is " + ch);
					chain.add(aacids.get(j));
					//chain.add(new AAcid(aacids.get(j).name(), aacids.get(j).hydrophobicity(), aacids.get(j).polarity(), aacidLocations[i]));
					//	System.out.println("at chain position " + chain.get(chain.size()-1).chainPosition() + " aminoacid " + chain.get(chain.size()-1).name());
					break;
				}
			}
		}
		return chain;
	}

	public static String pdbToString(String filename){
		String pdbString = "";

		//reading from atom doesnt work because it has missing sequence. have to read from seqres. atom part has information containing the structure.
		/*
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

		 */



		try {
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
					//System.out.println("a is " + a + " b is " + b);
					ch[a] = ch[a] + CodeToLetter(aStr);
					aNoS[a][b] = aNo;
					//System.out.println("at aminoacid #" + aNo + " there is " + aStr + " first letter " + CodeToLetter(aStr) + " line " + line);
				}
			}

			sc2.close();
			aacidLocations = new int[ch[a].length()];
			int c = 0;
			for (int i = 0; i < aNoS[a].length; i++){ //this will only be used if there is an issue with AminoAcids Not matching with the right locations
				if (aNoS[a][i] != 0){
					//System.out.print((aacidLocations[c] = aNoS[a][i]) + "-");
					aacidLocations[c] = aNoS[a][i];
					c++;
				}
			}


		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return pdbString;
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

	public static ArrayList<ChainAAcid> sort(ArrayList<ChainAAcid> sortedTendencies, String fasta, ArrayList<AAcid> chain, 
			ArrayList<Integer> hydrophobicRegions, ArrayList<Integer> hydrogenBonds, ArrayList<Integer> hydrophilicRegions, 
			ArrayList<Integer> noHydrogenBonds, ArrayList<AAcid> aacids ){

		double tendency1[] = new double[fasta.length()]; //hydrophobicity
		double tendency2[] = new double[fasta.length()]; //polarity
		double tendency3[] = new double[fasta.length()]; //hydrophilicity
		double tendency4[] = new double[fasta.length()]; //non-polarity

		for(int i = 0; i < fasta.length(); i++){
			if(hydrophobicRegions.contains(aacidLocations[i]) ){
				tendency1[i] = (5.7 - chain.get(i).hydrophobicity())/18.0; // normalize tendencies
				//System.out.println("Aminoacid " + aacidLocations[i] + " which is " + chain.get(i).name() + " should be hydrophobic");
			}
			if(hydrogenBonds.contains(aacidLocations[i])){
				tendency2[i] = (52.0 - chain.get(i).polarity())/52.0;
				//System.out.println("Aminoacid " + aacidLocations[i] + " which is " + chain.get(i).name() + " should have hydrogen bonds");
			}
			if(hydrophilicRegions.contains(aacidLocations[i]) ){  //burasi bir muamma??
				tendency3[i] = (chain.get(i).hydrophobicity()-(-12.3))/18.0; // normalize tendencies
				//System.out.println("Aminoacid " + aacidLocations[i] + " which is " + chain.get(i).name() + " should be hydrophilic");
			}
			if(noHydrogenBonds.contains(aacidLocations[i])){
				tendency4[i] = (chain.get(i).polarity()-0)/52.0;
				//System.out.println("Aminoacid " + aacidLocations[i] + " which is " + chain.get(i).name() + " should have no hydrogen bonds");
			}


		}

		//there should be a formula that finds the balance between hydrophobicity and hydrogen bonding

		double tendency[] = new double[fasta.length()];
		for (int i = 0; i < fasta.length(); i++ ){
			tendency[i] = scalar1*tendency1[i] + scalar2*tendency2[i] + scalar1*tendency3[i] + scalar2*tendency4[i] ;
			//System.out.println("tendency of aminoacid " + aacidLocations[i] + " is " + tendency1[i] + " + " + tendency2[i] + " = " + tendency[i]);
		}

		//determine the greatest tendency
		int greatestIndex = 0;
		double greatest = 0.0;
		for (int i = 0; i < fasta.length(); i++ ){
			if(sortedTendencies.size()==0){
				sortedTendencies.add(new ChainAAcid(chain.get(i).name(), chain.get(i).hydrophobicity(), 
						aacids.get(i).polarity(), aacidLocations[i], tendency[i], tendency1[i], tendency2[i], tendency3[i], tendency4[i]));
				//System.out.println("added " + sortedTendencies.get(0).toString() );
				continue;
			}
			else{
				for(int j=0; j < sortedTendencies.size(); j++){
					if(sortedTendencies.get(j).tendency() > tendency[i]){
						sortedTendencies.add(j, new ChainAAcid(chain.get(i).name(), chain.get(i).hydrophobicity(), 
								chain.get(i).polarity(), aacidLocations[i], tendency[i], tendency1[i], tendency2[i], tendency3[i], tendency4[i]));
						//	System.out.println("added aacid to tendencies: " + sortedTendencies.get(j).name() + " with tendency " + sortedTendencies.get(j).tendency());
					//	System.out.println("added " + sortedTendencies.get(j).toString() );

						break;
					}
					else if(j == (sortedTendencies.size()-1)){
						sortedTendencies.add(new ChainAAcid(chain.get(i).name(), chain.get(i).hydrophobicity(), 
								chain.get(i).polarity(), aacidLocations[i], tendency[i], tendency1[i], tendency2[i], tendency3[i], tendency4[i]));
						//System.out.println("added aacid to tendencies: " + sortedTendencies.get(j+1).name() + " with tendency " + sortedTendencies.get(j+1).tendency());
						//System.out.println("added " + sortedTendencies.get(sortedTendencies.size()-1).toString() );
						break;

					}
					else{
						continue;
					}
				}
			}
			/*	System.out.println("sorted tendency length is now " + sortedTendencies.size());
			if(tendency[i] > greatest){
				greatest = tendency[i];
				greatestIndex = i; 
			}*/
		}
		//System.out.println("sorted tendencies are " + sortedTendencies.toString());
		//System.out.println("The change should be at " + aacidLocations[greatestIndex] + " the aminoacid: " + chain.get(greatestIndex).name() );

		return sortedTendencies;
	}

	public static String sendToCommandline(String path, String input){
		String output = "";
		try {
			Process tr = Runtime.getRuntime().exec(path  );
			Writer wr = new OutputStreamWriter( tr.getOutputStream() );
			BufferedReader rd = new BufferedReader( new InputStreamReader( tr.getInputStream() ) );
			if(input != null || input != ""){
				wr.write( input );
			}
			wr.flush();
			String s;
			while((s = rd.readLine()) != null)
			{
				output = output + s + "\n";			
			}

			int exitVal = tr.waitFor();
			System.out.println("Exited with error code "+exitVal);
			rd.close();
			wr.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return output;
	}

	public static String tendencyToString(ArrayList<ChainAAcid> sortedTendencies, char[] predicted){
		String output = "";
		String output2 = "";
		String outputCombined = "";
		if(possibilities >= sortedTendencies.size()){
			System.out.println("asked possibilities exceed the number of aminoacids");
			System.exit(1);
		}
		for(int i = 0; i < possibilities; i++){
			output = "[";
			output2 = "[";
			for (int j = 0; j <= i; j++){
				int pos = sortedTendencies.size()-1-j;
			//	System.out.println("change priority " + (j+1)  + " should be at " + sortedTendencies.get(pos));
				output = output + "(" + sortedTendencies.get(pos).chainPosition() + "," +  predicted[2*j] +  ")" ;  
				output2 = output2 + "(" + sortedTendencies.get(pos).chainPosition() + "," +  predicted[2*j+1] +  ")" ;
			}
			outputCombined = outputCombined + output + "]\n";
			outputCombined = outputCombined + output2 + "]\n";
		}
		System.out.println("output =\n" + outputCombined);
		return outputCombined;

	}

	public static char[] predictChange(ArrayList<ChainAAcid> sortedTendencies, char[] sortedHydrophobicity, char[] sortedPolarity){
		char[] predicted = new char[possibilities*2];
		for(int i = 0; i < possibilities; i++){
			ChainAAcid currentAAcid = sortedTendencies.get(sortedTendencies.size()-i-1); 
			boolean hydrophobic = (currentAAcid.hydrophobicTendency() >= currentAAcid.polarTendency()) && 
					(currentAAcid.hydrophobicTendency() >= currentAAcid.nonPolarTendency()) &&
					(currentAAcid.hydrophobicTendency() >= currentAAcid.hydrophilicTendency());

			boolean hydrophilic = (currentAAcid.hydrophilicTendency() >= currentAAcid.polarTendency()) && 
					(currentAAcid.hydrophilicTendency() >= currentAAcid.nonPolarTendency()) &&
					(currentAAcid.hydrophilicTendency() >= currentAAcid.hydrophobicTendency());
			
			boolean polar = (currentAAcid.polarTendency() >= currentAAcid.hydrophobicTendency()) && 
					(currentAAcid.polarTendency() >= currentAAcid.nonPolarTendency()) &&
					(currentAAcid.polarTendency() >= currentAAcid.hydrophilicTendency());
			
			boolean nonpolar = (currentAAcid.nonPolarTendency() >= currentAAcid.hydrophobicTendency()) && 
					(currentAAcid.nonPolarTendency() >= currentAAcid.polarTendency()) &&
					(currentAAcid.nonPolarTendency() >= currentAAcid.hydrophilicTendency());
			
					//	System.out.println("hydrophobic tendency at this point is " + currentAAcid.hydrophobicTendency());
					//	System.out.println("polar tendency at this point is " + currentAAcid.polarTendency());
					if(hydrophobic){  //the change should be for a more hydrophobic aminoacid
						int length = sortedHydrophobicity.length;
						//first, add the most hydropobic aminoacid
						predicted[2*i] = sortedHydrophobicity[length-1];
						//	System.out.println("first add the highest = " + 2*i + " at " + predicted[2*i]);
						//then, find the one in the middle and add it
						for(int j = 0; j < length; j++){
							if(currentAAcid.name() - sortedHydrophobicity[j] == 0){
								//this is where our aacid stands
								predicted[(2*i+1)] = sortedHydrophobicity[(j+ (length-j)/2)];
								//	System.out.println("our aminoacid is at " + j +" second add the middle = " + (2*i+1) + " at " + predicted[2*i+1] + " hphobicity " + currentAAcid.hydrophobicTendency());
								break;
							}
						}

					}
					else if(hydrophilic){
						int length = sortedHydrophobicity.length;
						predicted[2*i] = sortedHydrophobicity[0];
						for(int j = 0; j < length; j++){
							if(currentAAcid.name() - sortedHydrophobicity[j] == 0){
								predicted[(2*i+1)] = sortedHydrophobicity[j/2];
								break;
							}
						}
					}
					else if(nonpolar){
						int length = sortedPolarity.length;
						predicted[2*i] = sortedPolarity[0];
						for(int j = 0; j < length; j++){
							if(currentAAcid.name() - sortedPolarity[j] == 0){
								predicted[(2*i+1)] = sortedPolarity[j/2];
								break;
							}
						}						
					}
					else{ //polar
						int length = sortedPolarity.length;
						//first, add the most hydropobic aminoacid
						predicted[2*i] = sortedPolarity[length-1];
						//	System.out.println("first add the highest = " + 2*i + " at " + predicted[2*i]+ " polarity " + currentAAcid.polarTendency());
						//then, find the one in the middle and add it
						for(int j = 0; j < length; j++){
							if(currentAAcid.name() - sortedPolarity[j] == 0){
								//this is where our aacid stands
								predicted[(2*i+1)] = sortedPolarity[j+(length-j)/2];
								//	System.out.println("our aminoacid is at " + j + " second add the middle = " + (2*i+1) + " at " + predicted[2*i+1] + " at " + (j+ (length-j)/2));
								break;
							}
						}

					}
		}

		return predicted;
	}

	public static char[] sortHydrophobics(ArrayList<AAcid> aacids){
		ArrayList<AAcid> sorted = new ArrayList<AAcid>();
		for(int i = 0; i < aacids.size(); i++){
			if(sorted.size()==0){
				sorted.add(aacids.get(i));
			}
			else{
				for(int j=0; j < sorted.size(); j++){
					if(sorted.get(j).hydrophobicity() > aacids.get(i).hydrophobicity()){  //found the first aacid that has greater hydrophobicity than our aacid.
						sorted.add(j, aacids.get(i));
						break;
					}
					else if(j == (sorted.size()-1)){ //this is the aacid with greatest hydrophobicity so far
						sorted.add(aacids.get(i));
						break;

					}
					else{
						continue;
					}
				}
			}
		}



		char[] sortedhydrophobics = new char[sorted.size()];
		for(int i = 0; i < sorted.size(); i++){
			sortedhydrophobics[i] = sorted.get(i).name();
			//	System.out.println("sorted is "  + i + "- " +sorted.get(i).name() + " at hydrophobicity " + sorted.get(i).hydrophobicity());
		}

		return sortedhydrophobics;
	}

	public static char[] sortPolarity(ArrayList<AAcid> aacids){
		ArrayList<AAcid> sorted = new ArrayList<AAcid>();
		for(int i = 0; i < aacids.size(); i++){
			if(sorted.size()==0){
				sorted.add(aacids.get(i));
			}
			else{
				for(int j=0; j < sorted.size(); j++){
					if(sorted.get(j).polarity() > aacids.get(i).polarity()){  //found the first aacid that has greater polarity than our aacid.
						sorted.add(j, aacids.get(i));
						break;
					}
					else if(j == (sorted.size()-1)){ //this is the aacid with greatest polarity so far
						sorted.add(aacids.get(i));
						break;

					}
					else{
						continue;
					}
				}
			}
		}



		char[] sortedpolarity = new char[sorted.size()];
		for(int i = 0; i < sorted.size(); i++){
			sortedpolarity[i] = sorted.get(i).name();
			//	System.out.println("sorted is " + i + "- " +sorted.get(i).name() + " at polarity " + sorted.get(i).polarity());
		}

		return sortedpolarity;
	}

}


