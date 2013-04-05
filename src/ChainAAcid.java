
public class ChainAAcid extends AAcid {
	int position;
	double tendency;
	
	public ChainAAcid(char aacidn, double hphobicity, double htendency,
			int position, double tendency) {
		super(aacidn, hphobicity, htendency);
		this.position = position;
		this.tendency = tendency;
	}
	
	public int chainPosition(){
		return position;
	}
	
	public double tendency(){
		return tendency;
	}
	
	public String toString(){
		return "AAcid " + name() + " at position " + chainPosition() + " tendency " + tendency();
	}

}
