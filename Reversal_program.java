package program;

import java.io.*;

public class SortByRevfinal {
	static int MAX_GENEGAP = 3;
	static int MIN_VALUE = -10000;
	static int MAX_VALUE = 10000;

	static int NUMSpecies = 0;
	static int NUMGenes = 0;
	int[][] dataMatrix = null;
	BreakP[][] BPgenes = null;
	int[][] ADJgenes = null;
	String[] StrainLabels = null;

	class BreakP {
		int gene = 0;
		int count = 1;
		BreakP next = null;

		BreakP(int g, BreakP n) {
			gene = g;
			next = n;
		}
	}

	public void readData(String name) {
		// Scan the file and decide the matrix size first.
		try {
			BufferedReader br = new BufferedReader(new FileReader(name));
			String line = br.readLine();
			NUMSpecies = line.split("\t").length;
			NUMGenes = 1;
			while ((line = br.readLine()) != null) {
				String[] words = line.split("\t");
				for (int i = 0; i < words.length; i++) {
					if (words[i].equals(""))
						continue;
					int I = Integer.parseInt(words[i]);
					if (I > NUMGenes)
						NUMGenes = I;
				}
			}
			br.close();
			System.out.println("NUMgenes: " + NUMGenes);
			System.out.println("NUMSpecies: " + NUMSpecies);
			// Initialization
			dataMatrix = new int[NUMGenes][NUMSpecies];
			BPgenes = new BreakP[2][NUMGenes];
			ADJgenes = new int[2][NUMGenes];
			StrainLabels = new String[NUMSpecies];
			for (int i = 0; i < NUMSpecies; i++)
				StrainLabels[i] = new String(new byte[] { (byte) (65 + i) }, "US-ASCII");
			// Reopen the file to read its contents.
			br = new BufferedReader(new FileReader(name));
			int row = 0;
			while ((line = br.readLine()) != null) {
				String[] words = line.split("\t");
				// We decrease the gene ID by 1 to fit Java / C structure.
				for (int i = 0; i < NUMSpecies; i++) {
					if (words.length <= i || words[i] == null || words[i].equals(""))
						dataMatrix[row][i] = MAX_VALUE;
					else {
						dataMatrix[row][i] = Integer.parseInt(words[i]) - 1;
					}
				}
				row++;
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	BreakP doesBPexist(int i, BreakP link) {
		while (link != null) {
			if (link.gene == i)
				return link;
			link = link.next;
		}
		return null;
	}

	String toBPString(BreakP link) {
		String s = "";
		while (link != null) {
			s += link.gene + "(" + link.count + ") ";
			link = link.next;
		}
		return s;
	}

	int getMaxNeighboringGene(BreakP link, int i) {
		// find the most likely neighboring number for gene i.
		int ret = -1;
		int c = 0;
		while (link != null) {
			if (c < link.count) {
				c = link.count;
				ret = link.gene;
			} else if (c == link.count) {
				// for tie breaking, choose the closer one
				if (Math.abs(ret - i) >= Math.abs(link.gene - i))
					ret = link.gene;
			}
			link = link.next;
		}
		return ret;
	}

	void createAdjacency() {
		// register neighboring gene ids
		for (int i = 0; i < NUMSpecies; i++) {
			for (int j = 0; j < NUMGenes; j++) {
				int cur = dataMatrix[j][i];
				if (cur == MAX_VALUE)
					break;
				int prev = (j == 0) ? MIN_VALUE : dataMatrix[j - 1][i];
				int next = (j == NUMGenes - 1) ? MAX_VALUE : dataMatrix[j + 1][i];
				if (prev > next) {
					int tmp = prev;
					prev = next;
					next = tmp;
				}
				BreakP p = doesBPexist(prev, BPgenes[0][cur]);
				if (p == null)
					BPgenes[0][cur] = new BreakP(prev, BPgenes[0][cur]);
				else
					p.count++;
				p = doesBPexist(next, BPgenes[1][cur]);
				if (p == null)
					BPgenes[1][cur] = new BreakP(next, BPgenes[1][cur]);
				else
					p.count++;
			}
		}
		// majority voting
		for (int i = 0; i < NUMGenes; i++) {
			ADJgenes[0][i] = getMaxNeighboringGene(BPgenes[0][i], i);
			ADJgenes[1][i] = getMaxNeighboringGene(BPgenes[1][i], i);
		}
		// flip neighbors in ADJgenes for reversal positions
		int prev = 0;
		int pos = 0;
		for (int i = 0; i < NUMGenes; i++) {
			int next = ADJgenes[1][pos];
			if (next == prev) {
				next = ADJgenes[0][pos];
				int tmp = ADJgenes[0][pos];
				ADJgenes[0][pos] = ADJgenes[1][pos];
				ADJgenes[1][pos] = tmp;
			}
			prev = pos;
			pos = next;
			if (pos == MAX_VALUE)
				break;
		}
	}

	void reorder(boolean verbose) {
		// This function reorders all gene numbers.
		int[] trans = new int[NUMGenes];
		int newPos = 0;
		for (int i = 0; i < NUMGenes; i++) {
			int next = ADJgenes[1][newPos];
			trans[newPos] = i;
			newPos = next;
			if (next == MAX_VALUE) {
				newPos = i;
				break;
			}
		}
		System.out.println("Consensus genes (max number): " + (newPos+1));

		if (verbose) {
			for (int i = 0; i < NUMGenes; i++)
				System.out.println(i + " -> " + trans[i] + "("
						+ toBPString(BPgenes[0][i]) + ": " + toBPString(BPgenes[1][i])
						+ ")");
		}

		for (int i = 0; i < NUMSpecies; i++) {
			for (int j = 0; j < NUMGenes; j++) {
				if (dataMatrix[j][i] != MAX_VALUE)
					dataMatrix[j][i] = trans[dataMatrix[j][i]];
			}
		}
	}

	void printGenes() {
		String s = "";
		for (int j = 0; j < NUMGenes; j++) {
			for (int i = 0; i < NUMSpecies; i++) {
				s += " " + dataMatrix[j][i];
			}
			s += "\n";
		}
		System.out.println(s);
	}

	int[] getBreakpoints(int species) {
		int[] tmp = new int[NUMGenes];
		int c = 0;
		for (int i = 0; i < NUMGenes - 2; i++) {
			if (dataMatrix[i + 1][species] == MAX_VALUE)
				break;
			if (dataMatrix[i + 1][species] == 0)
				continue;
			if (Math.abs(dataMatrix[i][species] - dataMatrix[i + 1][species]) > MAX_GENEGAP) {
				// single deletion or insertion is NOT a breakpoint
				tmp[c++] = i;
				i++;
			}
		}
		int[] ret = new int[c];
		for (int i = 0; i < c; i++) {
			ret[i] = tmp[i];
		}
		//System.out.println(ret.length);//added
		return ret;
	}

	int[] getGenes(int[] bp, int species) {
		int[] ret = new int[bp.length];
		for (int i = 0; i < bp.length; i++)
			ret[i] = dataMatrix[bp[i]][species];
		return ret;
	}

	String getStrainLabel(int species) {
		return StrainLabels[species];
	}

	void printBreakpoints(int[] bp, int species, boolean verbose) {
		// we assume the first position (gene 1) is not a breakpoint.
		StringBuffer sb = new StringBuffer();
		sb.append("strain ");
		sb.append(getStrainLabel(species));
		sb.append(" ... ");
		for (int i = 0; i < bp.length; i++) {
			int j = bp[i];
			if (!verbose) {
				sb.append(' ');
				sb.append(dataMatrix[j][species]);
			} else {
				if (j > 0) {
					sb.append(dataMatrix[j - 1][species]);
					sb.append(' ');
				}
				sb.append(dataMatrix[j][species]);
				if (dataMatrix[j + 1][species] != MAX_VALUE) {
					sb.append(" | ");
					sb.append(dataMatrix[j + 1][species]);
					sb.append(' ');
					sb.append(dataMatrix[j + 2][species]);
					sb.append(" ... ");
				}
			}
		}
		System.out.println(sb.toString());
	}

	boolean isSameBreakpoints(int[] x, int[] y) {
		if (x.length != y.length)
			return false;
		for (int i = 0; i < x.length; i++) {
			if (x[i] != y[i])
				return false;
		}
		return true;
	}

	void reduceSpecies() {
		// scan all breakpoints and reduces the species number
		for (int i = 0; i < NUMSpecies - 1; i++) {
			int[] breaki = getGenes(getBreakpoints(i), i);
			for (int j = i + 1; j < NUMSpecies; j++) {
				int[] breakj = getGenes(getBreakpoints(j), j);
				if (isSameBreakpoints(breaki, breakj)) {
					System.out.println("merge " + getStrainLabel(i) + " and "
							+ getStrainLabel(j));
					// update the strainlabel
					StrainLabels[i] += StrainLabels[j];
					//System.out.println(StrainLabels[i]);   //added  this
					// remove the column J
					if (j != NUMSpecies - 1) {
						//System.out.println(j);   //added this
						StrainLabels[j] = StrainLabels[NUMSpecies - 1];
						//System.out.println(StrainLabels[j]);   //added this
						for (int k = 0; k < NUMGenes; k++) {
							dataMatrix[k][j] = dataMatrix[k][NUMSpecies - 1];
						}
						j--;
					}
					NUMSpecies--;
				}
			}
		}
	}

	boolean removeRareReversals() {
		boolean ret = false;
		// collect breakpoints
		int[] count = new int[NUMGenes];
		int[] species = new int[NUMGenes];
		int[] datapos = new int[NUMGenes];
		for (int i = 0; i < NUMSpecies; i++) {
			int[] bp = getBreakpoints(i);
			for (int j = 0; j < bp.length; j++) {
				int gene = dataMatrix[bp[j]][i];
				count[gene]++;
				species[gene] = i;
				datapos[gene] = bp[j];
				//System.out.println(bp[j]);
				
			}
		}
		
		// locate strain-specific points
		for (int i = 0; i < NUMGenes - 1; i++) {
			if ((count[i] != 1) || count[i + 1] != 1)     
				continue;
			if (species[i] != species[i + 1])
				continue;
			// Check adjacency.
			if (Math.abs(dataMatrix[datapos[i]+1][species[i]]
					- dataMatrix[datapos[i + 1]+1][species[i]]) > MAX_GENEGAP)     //changed here removed +1 outside datapos[i]+1 and datapos[i+1]+1 
				continue;
			// Reverse genes between datapos[i]+1 and datapos[i+1]
			System.out.println(getStrainLabel(species[i]) + " Reversal " + i + " "
					+ (i + 1));
			ret = true;
			//System.out.println(dataMatrix.length);   //added this
			for (int j = 0; j < (Math.abs(datapos[i + 1] - datapos[i])) / 2; j++) {
				/*//System.out.println(datapos[i] + 1 + j);
				if ((datapos[i] + 1 + j)<dataMatrix.length){
					if(datapos[i]<datapos[i+1]){   //added this
				int tmp = dataMatrix[datapos[i] + 1 + j][species[i]];
				dataMatrix[datapos[i] + 1 + j][species[i]] = dataMatrix[datapos[i + 1]- j][species[i]];  
				//System.out.println(dataMatrix[datapos[i] + 1 + j][species[i]]); //added this line
				dataMatrix[datapos[i + 1] - j][species[i]] = tmp;
				//System.out.println(dataMatrix[datapos[i + 1] - j][species[i]]); //added this line
					}else{      //added this else
						//System.out.println(j);
						int tmp =dataMatrix[datapos[i+ 1]+1+j][species[i]];
						//System.out.println(tmp);
						dataMatrix[datapos[i+1] +1 + j][species[i]]=dataMatrix[datapos[i] - j][species[i]];
						//System.out.println(dataMatrix[datapos[i] - j][species[i]]);
						dataMatrix[datapos[i] - j][species[i]]= tmp;
						
						
					}
						
			}
			}*/
				if(datapos[i]<datapos[i+1]){
					if ((datapos[i] + 1 + j)<dataMatrix.length){
							int tmp = dataMatrix[datapos[i] + 1 + j][species[i]];
							dataMatrix[datapos[i] + 1 + j][species[i]] = dataMatrix[datapos[i + 1]- j][species[i]];  
							//System.out.println(dataMatrix[datapos[i] + 1 + j][species[i]]); //added this line
							dataMatrix[datapos[i + 1] - j][species[i]] = tmp;
							//System.out.println(dataMatrix[datapos[i + 1] - j][species[i]]); //added this line
					}
				}else{
						if((datapos[i+1]+1+j)<dataMatrix.length){
							int tmp =dataMatrix[datapos[i+ 1]+1+j][species[i]];
							//System.out.println(tmp);
							dataMatrix[datapos[i+1] +1 + j][species[i]]=dataMatrix[datapos[i] - j][species[i]];
							//System.out.println(dataMatrix[datapos[i] - j][species[i]]);
							dataMatrix[datapos[i] - j][species[i]]= tmp;
						}
					}
		}
		}
		return ret;
	}

	void removeReversal() {
		//int c1=0;  //added
		boolean f = true;
		while (f == true) {
			for (int i = 0; i < NUMSpecies; i++)
				printBreakpoints(getBreakpoints(i), i, true);
			reduceSpecies();	
			f = removeRareReversals();
			
			/*if(f==true){    //added
				c1=c1+1;
				if (c1==2){
					for (int i = 0; i <NUMSpecies ; i++) {
						for(int j=0;j<NUMGenes;j++){
							System.out.println(dataMatrix[j][i]);
						}
					}
			}
			}//added
*/		}
		
		for (int i = 0; i < NUMSpecies; i++)
			printBreakpoints(getBreakpoints(i), i, true);
	}
	
	void printConsensus(){
		int pos = 0;
		for (int i = 0; i < NUMGenes; i++){
			int next = ADJgenes[1][pos];
			System.out.println((i+1) + ":\t" + (pos +1));
			pos = next;		
			if (pos == MAX_VALUE)
				break;
		}
			
	}
	
	void writeBreakpoints(String filename){    //added this
		try{
			BufferedWriter bw= new BufferedWriter(new FileWriter(filename));
			for (int i = 0; i < NUMSpecies; i++) {
				int[] bp = getBreakpoints(i);
				bw.write(getStrainLabel(i)+":");
				if(bp.length!=0){
					for (int j = 0; j < bp.length; j++) {
						bw.write(dataMatrix[bp[j]][i]+ " ");
					}
				}
				else{
						bw.write('0');
				}
						
				
				bw.newLine();
			}
			bw.close();
		}
		catch (IOException e){}
	}        //till here
	
	
	void writeOrder(String filename){    //added this
		PrintWriter out = null;
		try{
			out = new PrintWriter (new FileWriter( filename ));
			for (int i = 0; i < NUMSpecies; i++) {
				
				for (int j = 0; j < NUMGenes; j++) {
//					System.out.println(dataMatrix[j][i]);
					out.print(dataMatrix[j][i]+ " ");
				}
				out.println();
			}
			out.close();
		}
		catch (IOException e){}
	}        //till here

	public static void main(String[] args) {
		SortByRevfinal test = new SortByRevfinal();
        test.readData("all_strains.txt");
		test.createAdjacency();
		//test.printConsensus();
		test.reorder(false);
//		test.writeOrder("Con_order.txt");
	/*for (int i = 23; i <NUMSpecies; i++) {
		for(int j=0;j<NUMGenes;j++){
			System.out.println(test.dataMatrix[j][i]);
		}
	}*/
		test.removeReversal();
		//test.writeBreakpoints("outputbp1.txt"); //added this
	}
}
