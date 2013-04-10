
public class ChainAAcid extends AAcid {
	int position;
	double tendency, tendency1, tendency2;
	
	public ChainAAcid(char aacidn, double hphobicity, double htendency,
			int position, double tendency, double tendency1, double tendency2) {
		super(aacidn, hphobicity, htendency);
		this.position = position;
		this.tendency = tendency;
		this.tendency1 = tendency1;
		this.tendency2 = tendency2;
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
	
	public String toString(){
		return "AAcid " + name() + " at position " + chainPosition() + " tendency " + tendency();
	}

}
