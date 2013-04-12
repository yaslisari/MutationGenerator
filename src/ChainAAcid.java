
public class ChainAAcid extends AAcid {
	int position;
	double tendency, tendency1, tendency2 ,tendency3, tendency4 ;
	
	public ChainAAcid(char aacidn, double hphobicity, double htendency,
			int position, double tendency, double tendency1, double tendency2, double tendency3, double tendency4) {
		super(aacidn, hphobicity, htendency);
		this.position = position;
		this.tendency = tendency;
		this.tendency1 = tendency1;
		this.tendency2 = tendency2;
		this.tendency3 = tendency3;
		this.tendency4 = tendency4;
	}
	
	public int chainPosition(){
		return position;
	}
	
	public double tendency(){
		return tendency;
	}
	
	public double hydrophobicTendency(){
		return tendency1;
	}
	
	public double polarTendency(){
		return tendency2;
	}
	
	public double hydrophilicTendency(){
		return tendency3;
	}
	
	public double nonPolarTendency(){
		return tendency4;
	}
	
	public String toString(){
		return "AAcid " + name() + " at position " + chainPosition() + " tendency " + tendency();
	}

}
