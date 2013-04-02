
public class AAcid {
	char acid;
	double hphob, polar;
	
	public AAcid(char aacidn, double hphobicity, double htendency ){
		acid = aacidn;
		hphob = hphobicity;
		polar = htendency;
	}
	
	public char name(){
		return acid;
	}
	
	public double hydrophobicity(){
		return hphob;
	}
	
	public double polarity(){
		return polar;
	}

}
