package edu.rice.cs.bioinfo.programs.phylonet.algos.iterHeuristic;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.start.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.TemporalConstraints;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.GLASSInference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeWithBranchLengthProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.start.distance.JCDistance;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import sun.nio.ch.Net;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by doriswang on 3/26/17.
 */
public class InferOperator {

    protected static String _msdir = "/Users/doriswang/msdir/ms";
    protected static String _ASTRALdir = "/Users/doriswang/PhyloNet/tools/Astral/";
    protected static String _RAxMLdir = "/Users/doriswang/PhyloNet/tools/standard-RAxML-master/";
    private static String TREE_DIR = "/Users/doriswang/PhyloNet/Data/IIG/tree/";
    private static String SEQ_DIR = "/Users/doriswang/PhyloNet/Data/IIG/seq/";
    private static String _fastTree = "/Users/doriswang/PhyloNet/tools/FT/";
    private static int iteration;

    InferOperator(int it) {
        this.iteration = it;
    }

    //TODO
    public static Map<String, String> readAlignmentFromNexusFile(String inFile, int len) {
        Map<String, String> aln = new HashMap<String, String>();

        return aln;
    }

    //TODO
    public static boolean writeAlignmentNexus(Map<String, String> aln, String outFile, int firstSite, int len) {
        boolean success = false;

        return success;
    }
    public static void main(String[] args) throws  IOException, InterruptedException, ParseException{

        InferOperator ifo = new InferOperator(1000);
        double[] ll = new double[50];
        String[] trees = new String[50];
        ifo.getInitTree(trees, ll, 50);
        //RawDataProcessor r = new RawDataProcessor();
        //String t = r.addTaxonName("(((19:0.00100000050002909,20:0.00100000050002909)I16:1.152946639631439,(16:0.4096753842384416,(18:0.00100000050002909,17:0.00100000050002909)I15:0.4086753837384125)I14:0.15927194422111118)I13:0.845163387858571,(((9:0.4456858217614752,(11:0.00100000050002909,10:0.00100000050002909)I11:0.4446858212614461)I10:0.2592441371312471,((15:0.00100000050002909,12:0.00100000050002909)I9:0.025201791973950474,(14:0.00100000050002909,13:0.00100000050002909)I8:0.12644601676689493)I7:0.07312226482896261)I6:0.4433325706515572,(3:0.2118922703896251,(2:0.20397143689818112,((1:0.00100000050002909,8:0.00100000050002909)I2:0.04122271582395445,((4:0.00100000050002909,5:0.00100000050002909)I17:0.08989500411654797,(6:0.00100000050002909,7:0.00100000050002909)I18:0.1380665477899139)I0:0.0)I1:0.06490488860823812)I3:0.007920833491443997)I4:1.0803345097601746)I5:0.3048617690139027)I12;");
        //System.out.println(t);
        //r.doubleMap();
//        String[] gts = ifo.getRAxMLStartT(38);
//        System.out.println(gts[37]);
//        int[] t1 = new int[]{665,525,420,260,260,469,796,296,220,262,256,194,447,849,798,522,267,274,471};
//        ifo.loadSYR(2, t1, 26, new ArrayList<Alignment>());


//        String[] taxaName = new String[taxaNum];
//        for (int i = 0; i < taxaNum; i++) {
//            line = new String[seqLength.length];
//            taxaName[i] = (String) br.readLine().trim().split(":")[0];
//            String temp = (String) br.readLine().trim().split(":")[1];
//            if (temp.length() < 8000)
//                System.err.println("Loading failed at : " + i);
//            int index = 0;
//            for (int j = 0; j < seqLength.length; j++) {
//                int len = seqLength[j];
//                if (j == seqLength[j] - 1)
//                    line[j] = temp.substring(index);
//                else {
//                    line[j] = temp.substring(index, index + len);
//                    index += len;
//                }
//            }
//            strArray.add(line);
//        }
//        for(int i = 0; i < seqLength.length; i++){
//            locus = new HashMap<String, String>();
//            fullLocus = new HashMap<String, String>();
//            for(int j = 0; j< taxaNum; j++){
//                locus.put(taxaName[j], strArray.get(i)[j]);
//                fullLocus.put(taxaName[j], strArray.get(i)[j]);
//            }
//            locus.remove("0");
//            Alignment fullAln = new Alignment(fullLocus);
//            Alignment aln = new Alignment(locus);
//            trueSeq.add(aln);
//            fullSeq.add(fullAln);
        }
//        String streeFile = "/Users/doriswang/PhyloNet/Data/17-taxon/004/ST1/1/Tree/trueST.txt";
//        BufferedReader stReader = new BufferedReader(new FileReader(streeFile));
//        String trueST = (String) stReader.readLine().trim();
//        String seqPath = "/Users/doriswang/PhyloNet/Data/17-taxon/004/ST1/1/Seq/";
//        List trueGTS = new ArrayList<Tree>();
//        List trueSeq = new ArrayList<Alignment>();
//
//        //TODO : 6.5 : count bp Number to decide if continue
//        for (int ln = 0; ln < 30 ; ln++) {
//            String tree = stReader.readLine().trim();
//            trueGTS.add(tree);
//        }
//        stReader.close();
//        for (int ln = 0; ln < 30; ln++) {
//            BufferedReader sReader = new BufferedReader(new FileReader(seqPath + "seq_" + ln + "_2000.nex"));
//            int j = 0;
//            while (j < 20) {
//                sReader.readLine();
//                j++;
//            }
//            String[] taxaName = new String[17];
//            String seqs = "";
//            String line = "";
//            Map<String, String> locus = new HashMap<String, String>();
//            for (int i = 0; i < 17; i++) {
//                line = sReader.readLine().trim();
//                String[] temp = line.split(" ");
//                taxaName[i] = temp[0];
//                if (temp[1].equals(""))
//                    seqs = temp[2];
//                else
//                    seqs = temp[1];
//                locus.put(taxaName[i], seqs.substring(0, 1000));
//            }
//            //locus.remove("O");
//            trueSeq.add(new Alignment(locus));
//
//            sReader.close();
//        }
//        ifo.initByRAxML(trueSeq,1000,17,30);
        //get16Aln(trueSeq, seqlens, 16,lociNum);


    public Map<String, Integer> loadSeqLens() throws IOException{
        Map<String, Integer> locusLen = new HashMap<String, Integer>();
        //<name, index>
        Map<String, String> locusName = new HashMap<String, String>();
        String inputFile = "/Users/doriswang/PhyloNet/Data/syr/syr011.txt";
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/Data/syr/locus_Mapping.txt"));
        BufferedWriter bwl = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/Data/syr/locus_Len.txt"));

        for(int i = 0 ; i< 63; i++) {
            br.readLine();
        }
        for(int i = 0 ; i< 19; i++) {
            String temp = br.readLine().trim();
            String t[] = temp.split(" ");
            String lName = t[1];
            String[] tempLen = t[3].split("-");
            int start = Integer.valueOf(tempLen[0]);
            int end = Integer.valueOf((tempLen[1].split(";")[0]));
            System.out.print(String.valueOf(end - start+1)  + ",");
            locusLen.put(lName, end - start+1);
            locusName.put(lName, String.valueOf(i));
            bw.write(i + ":" + lName + "\n");
            bwl.write(i + ":" +  String.valueOf(end - start+1) + "\n");
            bw.flush();
        }
        bw.close();
        for(int i = 0 ; i< 14; i++) {
            br.readLine();
        }

        //Read Tree

        for(int k = 0 ; k < 19; k++) {
            String path = "/Users/doriswang/PhyloNet/Data/SYRData/" + k + "/";
            //ifo.isExitsPath(path);
            bw = new BufferedWriter(new FileWriter(path + "startTree.txt"));

            String temp = br.readLine().trim();
            int sIndex = temp.indexOf("]");
            String stString = temp.substring(sIndex+1);
            int index1 = 1;
            int index0 = 0;
            List<Integer> indexes = new ArrayList<Integer>();
            for (int i = 0; i < stString.length(); i++) {
                if (stString.charAt(i) == ')') {
                    indexes.add(i);
                }
            }
            for (int i = indexes.size() - 2; i >0; i--) {
                index0 = indexes.get(i);
                for (int j = 1; j < stString.length()-index0; j++) {
                    index1 = index0 + j;
                    if (stString.charAt(index1) == ':') {
                        break;
                    }
                }
                String newST = stString.substring(0, index0+1) + "I" + Integer.toString(i) + stString.substring(index1 );
                stString = newST;
            }
            bw.write(stString + "\n");
            bw.flush();
            bw.close();
        }
        return locusLen;
    }

    public Alignment loadLocus(int seqNum, int seqLength, int taxaNum, String inputFile) throws IOException, InterruptedException {
        //"seq_" + tid + "_" + seqLen +".nex";scale + "/seq_" + i + "_" + len + ".nex"
        Map<String, String> locus = new HashMap<String, String>();
        inputFile = inputFile + "seq_" + seqNum + "_" + seqLength + ".nex";
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String line = "";
        int j = 1;
        while (j < 21) {
            br.readLine();
            j++;
        }
        for (int i = 0; i < taxaNum; i++) {
            line = br.readLine().trim();
            String[] temp = line.split((" "));
            if (temp.length > 2) {
                temp[1] = temp[2];
            }
            String lociName = temp[0];
            String seq = "";
            seq = temp[1].substring(0, seqLength);
            locus.put(lociName, seq);
            //locus.remove("O");
        }
        Alignment aln = new Alignment(locus);
        isExitsPath(_RAxMLdir + seqNum );
        BufferedWriter phy = new BufferedWriter(new FileWriter(_RAxMLdir + seqNum + "/dna.phy"));
        phy.write(taxaNum + " " + seqLength + '\n');
        phy.flush();
        List<String> names = aln.getTaxaNames();
        Map<String, String> thisAln = aln.getAlignment();
        for (int i = 0; i < taxaNum; i++) {
            String name = names.get(i);
            String seq = thisAln.get(name);
            phy.write(name + '\t' + seq + '\n');
        }
        phy.flush();
        phy.close();

        return aln;
    }


    public Alignment loadRandomLocus(int seqNum, int seqLength, int taxaNum, String inputFile, int offset) throws IOException, InterruptedException {
        //"seq_" + tid + "_" + seqLen +".nex";scale + "/seq_" + i + "_" + len + ".nex"
        Map<String, String> locus = new HashMap<String, String>();
        inputFile = inputFile + "seq_" + String.valueOf(seqNum+offset) + "_" + seqLength + ".nex";
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String line = "";
        int j = 1;
        while (j < 21) {
            br.readLine();
            j++;
        }
        for (int i = 0; i < taxaNum; i++) {
            line = br.readLine().trim();
            String[] temp = line.split((" "));
            if (temp.length > 2) {
                temp[1] = temp[2];
            }
            String lociName = temp[0];
            String seq = "";
            seq = temp[1].substring(0, seqLength);
            locus.put(lociName, seq);
            //locus.remove("O");
        }
        Alignment aln = new Alignment(locus);
        isExitsPath(_RAxMLdir + seqNum );
        BufferedWriter phy = new BufferedWriter(new FileWriter(_RAxMLdir + seqNum + "/dna.phy"));
        phy.write(taxaNum + " " + seqLength + '\n');
        phy.flush();
        List<String> names = aln.getTaxaNames();
        Map<String, String> thisAln = aln.getAlignment();
        for (int i = 0; i < taxaNum; i++) {
            String name = names.get(i);
            String seq = thisAln.get(name);
            phy.write(name + '\t' + seq + '\n');
        }
        phy.flush();
        phy.close();

        return aln;
    }


    //load seq data from raxml
    public List<Alignment> loadSYRByRX(int lociNum, int[] seqLength, int taxaNum, List<Alignment> trueSeq) throws IOException {
        //"seq_" + tid + "_" + seqLen +".nex";scale + "/seq_" + i + "_" + len + ".nex"
        Map<String, String> locus = new HashMap<String, String>(); //20 taxa except 1
        Map<String, String> fullLocus = new HashMap<String, String>(); // 21 taxa include i(OG)
        List<Alignment> fullalns = new ArrayList<Alignment>();
        for(int l = 0; l<lociNum; l++){
            String inputFile = "/Users/doriswang/PhyloNet/tools/standard-RAxML-master/" + l + "/dna.phy";
            BufferedReader br = new BufferedReader(new FileReader(inputFile));
            String line = "";
            String[] temp;
            line = br.readLine().trim();
            temp = line.split(" ");
            for (int i = 0; i < taxaNum; i++) {
                line = br.readLine().trim();
                temp = line.split((" "));
                if (temp.length > 2) {
                    temp[1] = temp[2];
                }
                String lociName = temp[0];
                //String seq = temp[1].substring(0, seqLength);
                locus.put(lociName, temp[1]);
                fullLocus.put(lociName, temp[1]);
                locus.remove("O");
            }
            br.close();
            Alignment fullAln = new Alignment(fullLocus);
            Alignment aln = new Alignment(locus);
            trueSeq.add(aln);
            fullalns.add(fullAln);
        }

        return fullalns;
    }


    // set true alignment
    // return full alignment
    public List<Alignment> loadSYR(int indivNum, int[] seqLength, int taxaNum, List<Alignment> trueSeq) throws IOException {
        //"seq_" + tid + "_" + seqLen +".nex";scale + "/seq_" + i + "_" + len + ".nex"
        Map<String, String> locus = new HashMap<String, String>();
        Map<String, String> fullLocus = new HashMap<String, String>();
        List<Alignment> fullSeq = new ArrayList<Alignment>();
        String inputFile = "/Users/doriswang/PhyloNet/Data/syr/" + indivNum+ "/seq.txt";
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String[] line = new String[seqLength.length];
        // List <taxon <seqs> >
        ArrayList<String[]> strArray = new ArrayList<String[]>();
        String[] taxaName = new String[taxaNum];
        for (int i = 0; i < taxaNum; i++) {
            line = new String[seqLength.length];
            String[] seqLines = br.readLine().trim().split(":");
            taxaName[i] = seqLines[0];
            String temp = seqLines[1];
            if (temp.length() < 8000)
                System.err.println("Loading failed at : " + i);
            int index = 0;
            for (int j = 0; j < seqLength.length; j++) {
                int len = seqLength[j];
                if (j == seqLength.length - 1)
                    line[j] = temp.substring(index);
                else {
                    line[j] = temp.substring(index, index + len);
                    index += len;
                }
            }
            strArray.add(line);
        }

        for(int i = 0; i < seqLength.length; i++){
            locus = new HashMap<String, String>();
            fullLocus = new HashMap<String, String>();
            for(int j = 0; j< taxaNum; j++){
                locus.put(taxaName[j], strArray.get(j)[i]);
                fullLocus.put(taxaName[j], strArray.get(j)[i]);
            }
            locus.remove("0");
            Alignment fullAln = new Alignment(fullLocus);
            Alignment aln = new Alignment(locus);
            trueSeq.add(aln);
            fullSeq.add(fullAln);
            BufferedWriter phy = new BufferedWriter(new FileWriter( _RAxMLdir + (i+19) + "/dna.phy" ));
            phy.write( taxaNum + " " + seqLength[i] + '\n');
            phy.flush();
            List<String> names = fullAln.getTaxaNames();
            Map<String, String> thisAln = fullAln.getAlignment();
            for (int j = 0; j < taxaNum; j++) {
                String name = names.get(j);
                String seq = thisAln.get(name);
                phy.write(name + '\t' + seq + '\n');
            }
            phy.flush();
            phy.close();
        }
        return fullSeq;
    }

    //Input: txt
    //Output: i/dna.phy


//  outgroup should be the child of root
    public String removeOutgroup(String tree) {
        //reroot
        // remove
        MutableTree thisTree = (MutableTree) Trees.readTree(tree.split(";")[0]);
        TMutableNode root = thisTree.getRoot();
        Iterator it = root.getChildren().iterator();
        TMutableNode otherChild = root;
        while (it.hasNext()) {
            TMutableNode child = (TMutableNode) it.next();
            if (child.getName().equals("O")) {
                if (it.hasNext())
                    otherChild = (TMutableNode) it.next();

                root.removeChild(child,false);
                thisTree.rerootTreeAtNode(otherChild);
                otherChild.removeChild(root,true);
                return thisTree.toString();
            }
            else
                otherChild = child;
        }
        return thisTree.toString();

    }

    //String seqFileName = "/Users/doriswang/PhyloNet/Data/IIG/seq/yeast_infer_network_from_multilocus_data_bayesian.nexus";
////t = 7 l = 106
    public List<Alignment> loadYData(String fileName, int taxaNum, int lociNum) throws IOException {
        List<Alignment> alns = new ArrayList<Alignment>();
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        int in = 0;
        //int[] lens = new int[106];
        while (in < 8) {
            br.readLine();
            in++;
        }
        int len = 0;
        for (int ln = 0; ln < lociNum; ln++) {
            Map<String, String> locus = new HashMap<String, String>();
            //List<Alignment> alns = new ArrayList<Alignment>();
            String line = "";
            String lociName = "";
            String seq = "";
            for (int i = 0; i < taxaNum; i++) {
                line = br.readLine().trim();

                String[] temp = line.split(" ");
                int index = 0;

                if (temp[index].length() > 1)
                    lociName = temp[0];
                else
                    lociName = temp[++index];
                if (temp[++index].length() > 1)
                    seq = temp[index];
                else
                    seq = temp[++index];
                locus.put(lociName, seq);
                //System.out.println("Locus_" + i +"'s length is " + seq.length());
                //len = seq.length();
            }
            //System.out.println("Locus_" + ln +"'s length is " + seq.length());
            System.out.print(" ," + seq.length());
            len += seq.length();
            //lens[lociNum] = seq.length();
            Alignment aln = new Alignment(locus);
            alns.add(aln);
            line = br.readLine();
            //if(line)
        }
        System.out.println("Ave length is " + len / 106);
        return alns;
    }


    //String seqFileName = "/Users/doriswang/PhyloNet/Data/IIG/seq/yeast_infer_network_from_multilocus_data_bayesian.nexus";
////t = 7 l = 106
    public List<Alignment> loadRData(String filePath, int repNum, int lociNum, List<Alignment> trueSeq, List<String> trueGTS) throws IOException {
        //#lociNum = 25
        // count bp Number to decide if continue
        for (int ln = 0; ln < lociNum; ln++) {
            String treeFile = filePath + "Rep" + repNum + ".gt" + ln + ".rose.tree.t";
            String seqFile = filePath + "Rep" + repNum + ".gt" + ln + ".rose.true.aln";
            BufferedReader tReader = new BufferedReader(new FileReader(treeFile));
            String tree = tReader.readLine().trim();
            trueGTS.add(tree);
            String[] lociName = new String[100];
            String[] seqs = new String[100];
            BufferedReader sReader = new BufferedReader(new FileReader(seqFile));
            String[] lines = sReader.readLine().trim().split(" ");
            //int tNum = Integer.getInteger(lines[0]);
            int tLength = Integer.valueOf(lines[1]);
            String line = "";
            for (int i = 0; i < 100; i++) {
                line = sReader.readLine().trim();
                String[] temp = line.split("    ");
                lociName[i] = temp[0];
                seqs[i] = temp[2].trim();
            }
            while (tLength != seqs[0].length()) {
                //SBF         GTGAATCAGGACGAAATTACTAGTAATACGCTAACAAGAATTCATGA--TATAGTGTCAA
                sReader.readLine();
                for (int i = 0; i < 100; i++) {
                    String l = sReader.readLine().trim();
                    seqs[i] += l;
                }
            }
            Map<String, String> locus = new HashMap<String, String>();
            for (int i = 0; i < 100; i++) {
                locus.put(lociName[i], seqs[i]);
            }
            trueSeq.add(new Alignment(locus));
        }
        //System.out.println("Ave length is " + len / 106);
        return trueSeq;
    }


//    //TODO: input: gtIndex: gtNum
//    public String loadRandomData (int stIndex, int seqlens, int lociNum, List<Alignment> trueSeq, List<String> trueGTS) throws IOException, ParseException, InterruptedException {
//        //#lociNum = 32
//        String streeFile = "/Users/doriswang/PhyloNet/Data/17-taxon/st/Rep" + stIndex +  "stree.txt";
//        BufferedReader stReader = new BufferedReader(new FileReader(streeFile));
//        String trueST = (String) stReader.readLine().trim();
//        String gtFile = "/Users/doriswang/PhyloNet/Data/17-taxon/32loci/Rep" + stIndex +  "gtrees.txt";
//        BufferedReader gtReader = new BufferedReader(new FileReader(streeFile));
//
//        for (int ln = 0; ln < lociNum ; ln++) {
//            String tree = stReader.readLine().trim();
//            trueGTS.add(tree);
//        }
//        stReader.close();
//
//        for (int ln = 0; ln < lociNum ; ln++) {
//            int tempLN = ln ;
//            BufferedReader sReader = new BufferedReader(new FileReader(seqPath + "seq_" + tempLN + "_2000.nex"));
//            int j = 0;
//            while (j < 20) {
//                sReader.readLine();
//                j++;
//            }
//            String[] taxaName = new String[17];
//            String seqs = "";
//            String line = "";
//            Map<String, String> locus = new HashMap<String, String>();
//            for (int i = 0; i < 17; i++) {
//                line = sReader.readLine().trim();
//                String[] temp = line.split(" ");
//                taxaName[i] = temp[0];
//                if (temp[1].equals(""))
//                    seqs = temp[2];
//                else
//                    seqs = temp[1];
//                locus.put(taxaName[i], seqs.substring(0, seqlens));
//            }
//            //locus.remove("O");
//            trueSeq.add(new Alignment(locus));
//            sReader.close();
//        }
//        int taxaNum = 17;
//        int seqLen = 1000;
//        for(int i = 0; i<lociNum; i++){
//            BufferedWriter phy = new BufferedWriter(new FileWriter(_RAxMLdir + i + "/dna.phy"));
//            phy.write( taxaNum + " " + seqLen + '\n');
//            phy.flush();
//            Alignment a = trueSeq.get(i);
//            List<String> names = a.getTaxaNames();
//            Map<String, String> thisAln = a.getAlignment();
//            for (int j = 0; j < taxaNum; j++) {
//                String name = names.get(j);
//                String seq = thisAln.get(name);
//                phy.write(name + '\t' + seq + '\n');
//            }
//            phy.flush();
//            phy.close();
//        }
//        get16Aln(trueSeq, seqlens, 16,lociNum);
//        return trueST;
//    }

    //String seqFileName = "/Users/doriswang/PhyloNet/Data/17-taxon/32loci/Rep467gtrees";
////t = 16 + 1 l = 32 length = 2000
    //TODO 7.18 : only fit for 17 taxon now!! need to be extened
    public String load17Data (String filePath, int seqlens, int lociNum, List<Alignment> trueSeq, List<String> trueGTS) throws IOException, ParseException, InterruptedException {
        //#lociNum = 32
        String streeFile = filePath + "Tree/trueST.txt";
        BufferedReader stReader = new BufferedReader(new FileReader(streeFile));
        String trueST = (String) stReader.readLine().trim();
        String seqPath = filePath + "Seq/";
//        for(int t = 0;t<10;t++){
//            stReader.readLine().trim();
//        }
        // 6.5 : count bp Number to decide if continue
        stReader = new BufferedReader(new FileReader(filePath + "Tree/trueGTS.txt"));
        for (int ln = 0; ln < lociNum ; ln++) {
            String tree = stReader.readLine().trim();
            trueGTS.add(tree);
        }
        stReader.close();
        for (int ln = 0; ln < lociNum ; ln++) {
            int tempLN = ln ;
            //TODO different file name
            BufferedReader sReader = new BufferedReader(new FileReader(seqPath + "seq_" + tempLN + "_1000.nex"));
            int j = 0;
            while (j < 20) {
                sReader.readLine();
                j++;
            }
            String[] taxaName = new String[17];
            String seqs = "";
            String line = "";
            Map<String, String> locus = new HashMap<String, String>();
            for (int i = 0; i < 17; i++) {
                line = sReader.readLine().trim();
                String[] temp = line.split(" ");
                taxaName[i] = temp[0];
                if (temp[1].equals(""))
                    seqs = temp[2];
                else
                    seqs = temp[1];
                locus.put(taxaName[i], seqs.substring(0, seqlens));
            }
            //locus.remove("O");
            trueSeq.add(new Alignment(locus));
            sReader.close();
        }
        int taxaNum = 17;
        int seqLen = seqlens;
        for(int i = 0; i<lociNum; i++){
            BufferedWriter phy = new BufferedWriter(new FileWriter(_RAxMLdir + i + "/dna.phy"));
            phy.write( taxaNum + " " + seqLen + '\n');
            phy.flush();
            Alignment a = trueSeq.get(i);
            List<String> names = a.getTaxaNames();
            Map<String, String> thisAln = a.getAlignment();
            for (int j = 0; j < taxaNum; j++) {
                String name = names.get(j);
                String seq = thisAln.get(name);
                phy.write(name + '\t' + seq + '\n');
            }
            phy.flush();
            phy.close();
        }
        get16Aln(trueSeq, seqlens, 16,lociNum);
        return trueST;
    }


    public static List<Tree> simulateGeneTrees(Network net,
                                               Map<String, List<String>> species2alleles, int numGT) {
        SimGTInNetworkByMS sim = new SimGTInNetworkByMS();
        return sim.generateGTs(Networks.readNetwork(net.toString()), species2alleles, numGT, _msdir);
    }

    //seq-gen file
    public List<Alignment> loadAlnSQ(String inputFile, int seqLength, int taxaNum, int loci) throws IOException {
        Map<String, String> locus = new HashMap<String, String>();
        List<Alignment> alns = new ArrayList<Alignment>();
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String line = "";
        int lineNum = 1;
        while (lineNum < 20) {
            br.readLine();
            lineNum++;
        }
        for (int j = 0; j < loci; j++) {
            if (j > 0) {
                int inter = 1;
                while (inter < 8) {
                    br.readLine();
                    inter++;
                }
            }
            for (int i = 0; i < taxaNum; i++) {
                line = br.readLine().trim();

                String[] temp = line.split(" ");
                if (temp.length > 2) {
                    temp[1] = temp[2];
                }
                String lociName = temp[0];
                String seq = "";
                seq = temp[1].substring(0, seqLength);
                locus.put(lociName, seq);
            }
            Alignment aln = new Alignment(locus);
            alns.add(aln);
        }

        return alns;
    }

    public static ArrayList<File> getListFiles(Object obj) {
        File directory = null;
        if (obj instanceof File) {
            directory = (File) obj;
        } else {
            directory = new File(obj.toString());
        }
        ArrayList<File> files = new ArrayList<File>();
        if (directory.isFile()) {
            files.add(directory);
            return files;
        } else if (directory.isDirectory()) {
            File[] fileArr = directory.listFiles();
            for (int i = 0; i < fileArr.length; i++) {
                File fileOne = fileArr[i];
                files.addAll(getListFiles(fileOne));
            }
        }
        return files;
    }

    public List<UltrametricTree> simGTByUPGMA(Alignment seq) {
        //BufferedReader stdout = new BufferedReader(new InputStreamReader(proc.getInputStream()));

        UPGMATree upgma = new UPGMATree(new JCDistance(seq.getAlignment()));
        UltrametricTree template = upgma.getUltrametricTree();
        System.out.println("Init gene tree： " + template.getTree().toNewick());
        List<UltrametricTree> geneTrees = new ArrayList<>();
        geneTrees.add(upgma.getUltrametricTree());
        return geneTrees;
    }

    public List<UltrametricTree> simGTSByUPGMA(List<Alignment> trueAlns) throws IOException {
        List<UltrametricTree> geneTrees = new ArrayList<UltrametricTree>();
        for (int i = 0; i < trueAlns.size(); i++) {
            geneTrees.addAll(simGTByUPGMA((Alignment) trueAlns.get(i)));
        }
        return geneTrees;
    }


    protected Tree getGLASSTree(List<UltrametricTree> gtList, Map<String, List<String>> species2alleles) {
        GLASSInference glassInference = new GLASSInference();
        Tree inferredST = glassInference.inferSpeciesTree(getTreeByUTree(gtList));
        return inferredST;
    }

    //
    public double getAverage(double[] array){
        int sum = 0;
        for(int i = 0;i < array.length;i++){
            sum += array[i];
        }
        return (double)(sum / array.length);
    }

    // P(gt|ST)
    //if not return -> add st height on all leaves
    public List<Double> getGTSLLBySTYF(List<Tree> geneTrees, Network<NetNodeInfo> currentST) {
        List<Double> gProST = new ArrayList<>();
        double[] result = new double[geneTrees.size()];
        Network<NetNodeInfo> ultraST = getUltraST(currentST.toString());
        GeneTreeWithBranchLengthProbabilityYF g = new GeneTreeWithBranchLengthProbabilityYF(ultraST, geneTrees, null);
        g.calculateGTDistribution(result);
        for (int i = 0; i < result.length; i++)
            gProST.add(Math.log(result[i]));

        return gProST;
    }

    //
    public double getStandardDevition(double[] array){
        double sum = 0;
        double ave = getAverage(array);
        for(int i = 0;i < array.length;i++){
            sum += Math.sqrt(((double)array[i] - ave) * (array[i] - ave));
        }
        return (sum / (array.length - 1));
    }
    public BniNetwork<NetNodeInfo> inferSTByGLASS(List<UltrametricTree> geneTrees, double halfTheta) throws IOException, ParseException {
        Tree speciesTree = getGLASSTree(geneTrees, null);
        Utils._ESTIMATE_POP_SIZE = false;
        Utils._CONST_POP_SIZE = true;
        Utils._POP_SIZE_MEAN = 0.036;
        Map<String, Double> constraints = TemporalConstraints.getTemporalConstraints(geneTrees, null, null);
        UltrametricNetwork stWithPara = new UltrametricNetwork(speciesTree.toNewick(), geneTrees, null);
        //stWithPara.initNetHeights(4*scale,constraints);
        //Networks.autoLabelNodes((Network)stWithPara.getNetwork());
        BniNetwork<NetNodeInfo> st = (BniNetwork) stWithPara.getNetwork();
        Networks.autoLabelNodes(st);
        return st;
    }

    public Tree inferSTByGLASS(List<Tree> geneTrees) throws IOException, ParseException {
        GLASSInference glassInference = new GLASSInference();
        Tree inferredST = glassInference.inferSpeciesTree(geneTrees);
        Utils._ESTIMATE_POP_SIZE = false;
        Utils._CONST_POP_SIZE = true;
        Utils._POP_SIZE_MEAN = 0.001;
        List<UltrametricTree> uList = getUTreeByTree(geneTrees);
        Map<String, Double> constraints = TemporalConstraints.getTemporalConstraints(uList, null, null);
        UltrametricNetwork stWithPara = new UltrametricNetwork(inferredST.toNewick(), uList, null);
        //stWithPara.initNetHeights(4*scale,constraints);
        //Networks.autoLabelNodes((Network)stWithPara.getNetwork());
        BniNetwork<NetNodeInfo> st = (BniNetwork) stWithPara.getNetwork();
        Networks.autoLabelNodes(st);
        Tree glassTree = Trees.readTree(st.toString());
        return glassTree;
    }

    //infer ST by GTS by ASTRAL
    public BniNetwork<NetNodeInfo> inferSTByAST(List<UltrametricTree> geneTrees, double popSize) throws IOException, ParseException {
        String command = "java -jar " + _ASTRALdir + "astral.4.11.1.jar -i /Users/doriswang/PhyloNet/tools/Astral/test_data/tempIn.tre -o /Users/doriswang/PhyloNet/tools/Astral/test_data/testOut.tre";
        String cmdFile = _ASTRALdir + "tempCMD.sh";
        BufferedWriter cmd = new BufferedWriter(new FileWriter(cmdFile));
        cmd.write(command + '\n');
        cmd.flush();
        cmd.close();
        String inputFile = _ASTRALdir + "test_data/tempIn.tre";
        BufferedWriter in = new BufferedWriter(new FileWriter(inputFile));
        for (int i = 0; i < geneTrees.size(); i++) {
            in.write(geneTrees.get(i).toString() + "\n");
        }
        in.flush();
        in.close();
        ProcessBuilder pb = new ProcessBuilder("/bin/bash", cmdFile);
        pb.redirectErrorStream(true);
        try {
            Process proc = pb.start();
            try {
                proc.waitFor();
            } catch (InterruptedException ex) {
                System.out.println(ex.getMessage());
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        String outputFile = _ASTRALdir + "test_data/testOut.tre";
        BufferedReader out = new BufferedReader(new FileReader(outputFile));
        String stString = out.readLine().trim();
        BniNetwork<NetNodeInfo> st1 = new BniNetwork<>();
        int index0 = 0;
        int index1 = 1;
        int count = 0;
        List<Integer> indexes = new ArrayList<Integer>();
        for (int i = 0; i < stString.length(); i++) {
            if (stString.charAt(index1) == ':') {
                indexes.add(i);
            }
        }
        for (int i = 1; i <= indexes.size(); i++) {
            index1 = indexes.get(indexes.size() - i);
            for (int j = 1; j < stString.length(); j++) {
                index0 = index1 - j;
                if (stString.charAt(index0) == ')') {
                    break;
                }
            }
            String newST = stString.substring(0, index0) + "I" + Integer.toString(i - 1) + stString.substring(index1 + 1);
            stString = newST;
        }
//        ArrayList<Double> paras = new ArrayList<Double>();
//        String topo = "";
//        String firstNode = "";
//        //String pattern1 = "^\\d+";
//        String pattern = "\\d+(\\.\\d+)?";
//        Pattern p = Pattern.compile(pattern);
//        String[] s1 = stString.split(":");
//        for(int i = 0; i<s1.length; i++){
//            Matcher m = p.matcher(s1[i]);
//            boolean hasStart = false;
//            int firstEnd = -1;
//            String temp = "";
//            while(m.find()){
//                int start = m.start();
//                int end = m.end();
//                if(start==0) {
//                    temp += s1[i].substring(end);
//                    hasStart = true;
//                    firstEnd = end;
//                    continue;
//                }
//                if(end==s1[i].length()) {
//                    if(hasStart)
//                        temp = s1[i].substring(firstEnd,start);
//                    else
//                        temp = s1[i].substring(0, start);
//                }
//            }
//            topo+=temp;
//        }
//
//        BniNetwork st = (BniNetwork) Networks.readNetwork(topo);
//        Networks.autoLabelNodes(st);
//        Utils._ESTIMATE_POP_SIZE = false;
//        Utils._CONST_POP_SIZE = true;
//        Utils._POP_SIZE_MEAN = popSize;
//        Map<String, Double> constraints = TemporalConstraints.getTemporalConstraints(geneTrees, null, null);
//        UltrametricNetwork stWithPara = new UltrametricNetwork(st.toString(), geneTrees, null);
//        //stWithPara.initNetHeights(4*scale,constraints);
//        //Networks.autoLabelNodes((Network)stWithPara.getNetwork());
//        BniNetwork<NetNodeInfo> st1 = (BniNetwork) stWithPara.getNetwork();
//        Networks.autoLabelNodes(st1);
//        return st1;

        return st1;
    }

    //./raxmlHPC -M -m GTRGAMMA -p 12345 -q simpleDNApartition.txt -s dna.phy -n T22
    public List<Tree> initByRAxML(List<Alignment> aln, int seqLen, int taxaNum, int lociNum) throws IOException, ParseException, InterruptedException {
        //writeSeq + partition     rm + raxml    read in
        List<Tree> gts = new ArrayList<Tree>();
        for (int i = 0; i < lociNum; i++) {
            System.out.println("RAXML No" + i + "Start...");
            isExitsPath(_RAxMLdir + i);
            String rm = "rm *.T" + i + '\n';
            //./raxmlHPC-PTHREADS -M -m GTRGAMMA -p 12345 -# 10 -s 0/dna.phy -n T0
            String raxml = _RAxMLdir + "raxmlHPC-PTHREADS -M -m GTRGAMMA -p 12345 -# 30 " + "-s " + i + "/dna.phy -n T" + i;
            String cmdFile = _RAxMLdir + "/" + i + "/tempCMD.sh";
            BufferedWriter cmd = new BufferedWriter(new FileWriter(cmdFile));
            cmd.write("cd " + _RAxMLdir.substring(0, _RAxMLdir.length() - 1) + '\n');
            cmd.write(rm);
            cmd.flush();
            cmd.write(raxml);
            cmd.flush();
            cmd.close();

            ProcessBuilder pb = new ProcessBuilder("/bin/bash", cmdFile);
            pb.redirectErrorStream(true);
            try {
                Process proc = pb.start();
                try {
                    proc.waitFor();
                } catch (InterruptedException ex) {
                    System.out.println(ex.getMessage());
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        for (int i = 0; i < lociNum; i++) {
            String outputFile = _RAxMLdir + "RAxML_bestTree.T" + i;
            BufferedReader gtReader = new BufferedReader(new FileReader(outputFile));
            String gtString = gtReader.readLine().trim();
            gts.add(Trees.readTree(gtString));
            gtReader.close();
        }
        return gts;
    }
//FastTree -gtr -nt < alignment.file > tree_file
    public List<Tree> initFasttree(List<Alignment> aln, int[] seqLen, int lociNum) throws IOException, ParseException, InterruptedException {
        //writeSeq + partition     rm + raxml    read in
        List<Tree> gts = new ArrayList<Tree>();
        for (int i = 0; i < lociNum; i++) {
            //System.out.println("No" + i + "Start...");
            BufferedWriter bw = new BufferedWriter(new FileWriter(_fastTree + i + ".phy"));
            Alignment temp = aln.get(i);
            Map<String, String> thisAln = temp.getAlignment();
            int taxaNum = temp.getTaxonSize();
            bw.write(taxaNum + " " + seqLen[i] + '\n');
            for (int j = 0; j < taxaNum; j++) {
                String name = temp.getTaxaNames().get(j);
                String seq = thisAln.get(name);
                bw.write(name + ' ' + seq + '\n');
            }
            bw.flush();
            bw.close();
            isExitsPath(_fastTree + i);
            String rm = "rm *.T" + i + '\n';
            //./raxmlHPC-PTHREADS -M -m GTRGAMMA -p 12345 -# 10 -s 0/dna.phy -n T0
            String command = _fastTree + "FastTree -gtr -nt < " + i + ".phy > " + i + ".T";
            String cmdFile = _fastTree + "FastTree" + i + ".sh";
            BufferedWriter cmd = new BufferedWriter(new FileWriter(cmdFile));
            cmd.write("cd " + _fastTree.substring(0, _fastTree.length() - 1) + '\n');
            cmd.write(rm);
            cmd.flush();
            cmd.write(command);
            cmd.flush();
            cmd.close();

            ProcessBuilder pb = new ProcessBuilder("/bin/bash", cmdFile);
            pb.redirectErrorStream(true);
            try {
                Process proc = pb.start();
                try {
                    proc.waitFor();
                } catch (InterruptedException ex) {
                    System.out.println(ex.getMessage());
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        for (int i = 0; i < lociNum; i++) {
            String outputFile = _fastTree + i + ".T";
            BufferedReader out = new BufferedReader(new FileReader(outputFile));
            //String gtString = gtReader.readLine().trim();
            String stString = out.readLine().trim();
            int index0 = 0;
            int index1 = 1;
            for(int j = 1; j< stString.length()-1;j++){
                if(stString.charAt(j)!=')'){
                    continue;
                }
                else if(stString.charAt(j+1)==':'){
                    continue;
                }
                else{
                    int offset = 1;
                    boolean find = false;
                    for(offset=1;offset<stString.length()-j;offset++){
                        if(stString.charAt(offset+j)==':') {

                            find = true;
                            break;
                        }
                        else
                            continue;
                    }
                    if(find){
                        stString = stString.substring(0,j+1) + stString.substring(offset+j);
                    }
                }
            }
            gts.add(Trees.readTree(stString));
            //System.out.println(stString);
            out.close();
        }
        return gts;
    }
    //TODO
    public String[] getRAxMLStartT(int lociNum) throws IOException, ParseException, InterruptedException {
        //List<Tree> gts = new ArrayList<Tree>();
        String[] gts = new String[lociNum];
        for (int i = 0; i < lociNum; i++) {
            String tree = _RAxMLdir + "/RAxML_Init/RAxML_bestTree.T" + i;
            String llFile = _RAxMLdir + "/RAxML_Init/RAxML_info.T" + i;
            BufferedReader gtReader = new BufferedReader(new FileReader(tree));
            String gtString = gtReader.readLine().trim();
            gtReader.close();
            BufferedReader llReader = new BufferedReader(new FileReader(llFile));
            double ll = 0.0;
            String llString = "";
            while(llString!=null){
                llString = llReader.readLine();
                if(llString.length()<20||llString.isEmpty())
                    continue;
                if(llString.substring(0,5).equals("Final")){
                    String[] result = llString.split(" ");
                    ll = Double.valueOf(result[6]);
                    break;
                }
            }
            gts[i] = gtString + " " + ll;
            llReader.close();
        }
        return gts;
    }


        //operator.getUpdateSh(lociNum, msSize);
    public void getUpdateSh(int lociNum, int msSize) throws IOException, ParseException, InterruptedException {
        //writeSeq + partition     rm + raxml    read in
        List<Tree> gts = new ArrayList<Tree>();
        for (int i = 0; i < lociNum; i++) {
            String seqPath = _RAxMLdir + i;
            for(int j = 0; j<msSize; j++){
                String tName = "L" + String.valueOf(i) + "T" + String.valueOf(j); // L_i T_j
                BufferedWriter phy = new BufferedWriter(new FileWriter(seqPath + "/updateL" + String.valueOf(i) + "T" + String.valueOf(j) + ".sh"));
                phy.write("cd /Users/doriswang/PhyloNet/tools/standard-RAxML-master" + '\n');
                phy.flush();
                phy.write("rm *." + tName + '\n');
                phy.flush();
                phy.write("/Users/doriswang/PhyloNet/tools/standard-RAxML-master/raxmlHPC-PTHREADS  -f e -t topo.T" + String.valueOf(j) + " -m GTRGAMMA -s " + String.valueOf(i) + "/dna.phy -n " + tName + '\n');
                phy.flush();
                phy.close();
            }

        }
        //return gts;
    }

    //write alignment with 16 taxa for RAxML
    public List<Tree> get16Aln(List<Alignment> aln, int seqLen, int taxaNum, int lociNum) throws IOException, ParseException, InterruptedException {
        //writeSeq + partition     rm + raxml    read in
        List<Tree> gts = new ArrayList<Tree>();
        for (int i = 0; i < lociNum; i++) {
            isExitsPath(_RAxMLdir + i);
            BufferedWriter phy = new BufferedWriter(new FileWriter(_RAxMLdir + i + "/dnaNoOG.phy"));
            phy.write(taxaNum + " " + seqLen + '\n');
            phy.flush();
            Alignment a = aln.get(i);
            List<String> names = a.getTaxaNames();
            Map<String, String> thisAln = a.getAlignment();
            for (int j = 0; j < taxaNum; j++) {
                String name = names.get(j);
                if (name == "O")
                    continue;
                String seq = thisAln.get(name);
                phy.write(name + '\t' + seq + '\n');
            }
            phy.flush();
            phy.close();
        }
        return gts;
    }

    //Undo: false -> *halftheta    true -> /halfTheta
    public static Network<NetNodeInfo> scaleSTNet(Network<NetNodeInfo> currentST, double ratio, boolean undo) {
        //Tree newT = new STITree<>();
        Iterator itNode = currentST.getTreeNodes().iterator();
        if (undo == false) {
            while (itNode.hasNext()) {
                NetNode n = (NetNode) itNode.next();
                if (n == (NetNode) currentST.getRoot())
                    continue;
                NetNode p = (NetNode) n.getParents().iterator().next();
                if (n.getParentDistance(p) != TMutableNode.NO_DISTANCE) {
                    n.setParentDistance(p, n.getParentDistance(p) * ratio); // height = height * ratio
                }
            }
        } else {
            while (itNode.hasNext()) {
                NetNode n = (NetNode) itNode.next();
                if (n == (NetNode) currentST.getRoot())
                    continue;
                NetNode p = (NetNode) n.getParents().iterator().next();
                if (n.getParentDistance(p) != TMutableNode.NO_DISTANCE) {
                    n.setParentDistance(p, n.getParentDistance(p) / ratio); // height = height / ratio
                }
            }
        }
        return currentST;
    }

    public boolean hasCheckedBL(BniNetNode node, Network st) {
        boolean hasChecked = true;
        while (!node.isLeaf()) {
            BniNetNode child = (BniNetNode) node.getChildren().iterator().next();
            if (child.getParentDistance(node) == Double.NEGATIVE_INFINITY)
                hasChecked = false;
            node = child;
        }
        return hasChecked;
    }

    //get distance(a,b) where common parent is cRoot on one gt
    public double getLocalDist(TNode left, TNode right, TNode cRoot, Tree gt) {
        double lDistance = left.getParentDistance();
        double rDistance = right.getParentDistance();

        TNode lParent = left.getParent();
        TNode rParent = right.getParent();
        while (lParent != cRoot) {
            lDistance += lParent.getParentDistance();
            lParent = lParent.getParent();
        }
        while (rParent != cRoot) {
            rDistance += rParent.getParentDistance();
            rParent = rParent.getParent();
        }
        return lDistance + rDistance;
    }
//    //p is child or child's ancestor
//    public boolean isAncestor(TNode parent, TNode child, Stack<TNode> s){
//        if(parent.equals(child))
//            return true;
//        else{
//            if(parent.isLeaf())
//                return false;
//            else {
//                Iterator it = parent.getChildren().iterator();
//                TNode temp1 = (TNode)it.next();
//                TNode temp2 = (TNode)it.next();
//                if(temp1.equals(child)||temp2.equals(child))
//                    return true;
//                else
//                    return isAncestor(temp1,child,s)|| isAncestor(temp2,child,s);
//            }
//        }
//    }


    //Input: true gene trees ; ms sampled trees for iteration
    //Output: the best distance IIG may get for each locus in each iteration
    public void getMSDist(List<Tree> trueTopos, List<Tree> msTrees, double[][] msTreeDist, int itNum){
        for(int i = 0; i<trueTopos.size();i++){
            Tree thistree = trueTopos.get(i);
            double minD = 1;
            for(int j = 0; j<msTrees.size();j++){
                Tree msT = msTrees.get(j);
                double d = getDistance(thistree,msT);
                if(d<minD)
                    minD = d;
            }
            msTreeDist[itNum][i] = minD;
        }
    }


    //get Lowest Common Parent on gt for computing distance(a,b)
    public TNode getLCP(TNode root, TNode p, TNode q) {
        //System.out.println("getLCP" + p.getName() + " -- " + q.getName());
        if (root == null)
            return null;
        if (root == p || root == q)
            return root;
        if (root.isLeaf())
            return null;
        else {
            Iterator it = root.getChildren().iterator();
            TNode l = getLCP((TNode) it.next(), p, q);
            TNode r = getLCP((TNode) it.next(), p, q);
            if (l != null && r != null) {
                return root;
            } else if (l == null && r == null) {
                return null;
            } else {
                return l == null ? r : l;
            }
        }
    }

    //getInternalHeight(sib, st, gts);
    public double getInternalHeight(NetNode thisNode, Network st, Map<String, Double> heightList) {
        double height = Double.NEGATIVE_INFINITY;
        Iterator children = thisNode.getChildren().iterator();
        while (children.hasNext()) {
            NetNode c = (NetNode) children.next();
            Double cHeight = 0.0;
            if(!heightList.containsKey(c.getName()))
                cHeight = getInternalHeight(c, st, heightList);
            if (height < cHeight + c.getParentDistance(thisNode)) {
                    height = heightList.get(c.getName()) + c.getParentDistance(thisNode);
                    heightList.put(thisNode.getName(), height);
                }
        }
        return height;

    }

    // return Map<#level, nodesList>
    public Map<Integer, List<String>> getLevels(Network st) {
        Map levelMap = new HashMap<Integer, List<String>>();
        LinkedList q = new LinkedList<NetNode>();
        q.offer(st.getRoot());
        int cur, last;
        int level = 1;
        NetNode current = new BniNetNode();
        // Map<Integer, String> map = new HashMap<Integer, String>();
        while (!q.isEmpty()) {
            cur = 0;
            last = q.size();
            while (cur < last) {
                current = (NetNode) q.poll();
                cur++;
                if (current.isLeaf()) {
                    if (levelMap.containsKey(level)) {
                        List nodes = (List) levelMap.get(level);
                        nodes.add(current.getName());
                    } else {
                        List<String> nodes = new ArrayList<String>();
                        nodes.add(current.getName());
                        levelMap.put(level, nodes);
                    }
                } else {
                    Iterator children = current.getChildren().iterator();
                    while (children.hasNext())
                        q.offer(children.next());
                }
            }
            level++;
        }
        return levelMap;
    }

    public Network<NetNodeInfo> getUltraST(String stString) {
        Network<NetNodeInfo> st = Networks.readNetwork(stString);
        List<NetNode> visit = new ArrayList<NetNode>();
        Map<NetNode,Double> heights = new HashMap<NetNode,Double>();
        Iterator leafIt = st.getLeaves().iterator();
        Double netHeight = 0.0;
        while(leafIt.hasNext()){
            NetNode temp = (NetNode)leafIt.next();
            if(visit.contains(temp))
                continue;
            else{
                visit.add(temp);
                double tempHeight = getLeafHeight(temp,st);
                if(tempHeight>netHeight)
                    netHeight = tempHeight;
                heights.put(temp,tempHeight);
                NetNode p = (NetNode) temp.getParents().iterator().next();
                NetNode sib = (NetNode) p.getChildren().iterator().next();
                if(sib.equals(temp))
                    sib = (NetNode) p.getChildren().iterator().next();
                if(sib.isLeaf()){
                    visit.add(sib);
                    heights.put(sib,tempHeight);
                }

            }
        }
        Iterator leaves = heights.keySet().iterator();
        while(leaves.hasNext()){
            NetNode thisNode = (NetNode)leaves.next();
            NetNode parent = (NetNode)thisNode.getParents().iterator().next();
            Double currentH = heights.get(thisNode);
            Double dist = thisNode.getParentDistance(parent);
            thisNode.setParentDistance(parent,dist + netHeight-currentH);

        }
        return st;
    }

    //    public double getLowestHeight(NetNode n, Network st){
//        double h = 0.0;
//        if(n.isLeaf())
//            return h;
//        Stack s = new Stack();
//        Iterator tempIt = n.getChildren().iterator();
//        s.add(n);
//        Map<String,Double> heights = new HashMap<String,Double>();
//        while(!s.isEmpty()){
//            NetNode thisNode = (NetNode) s.pop();
//            if(thisNode.isLeaf()){
//                heights.add
//            }
//        }
//        return h;
//    }
    public boolean isExitsPath(String path) throws InterruptedException {
        File file = new File(path.toString());
        if (!file.exists()) {
            file.mkdirs();
            System.out.println("build dir：" + path.toString());
            //Thread.sleep(1500);
        }
        File file1 = new File(path.toString());
        if (!file1.exists()) {
            return true;
        } else {
            return false;
        }
    }

    public void updateMSTopos(List<Tree> topos) throws IOException {

            String seqPath = _RAxMLdir;
            for(int j = 0; j<topos.size(); j++){
                BufferedWriter topow = new BufferedWriter(new FileWriter(seqPath + "topo.T" + String.valueOf(j)));
                topow.write(topos.get(j).toString());
                topow.flush();
                topow.close();
            }

    }


    public void runUpdateShell(int msSize, int lociNum) throws IOException {

        //String seqPath = _RAxMLdir;
        for (int i = 0; i < lociNum; i++) {
            String seqPath = _RAxMLdir + i + "/";
            for(int j = 0; j<msSize; j++){
                //updateL0T1.sh
                String cmdFile = _RAxMLdir + "/" + i + "/updateL" + String.valueOf(i) + "T" + String.valueOf(j) + ".sh";
                ProcessBuilder pb = new ProcessBuilder("/bin/bash", cmdFile);
                pb.redirectErrorStream(true);
                try {
                    Process proc = pb.start();
                    try {
                        proc.waitFor();
                    } catch (InterruptedException ex) {
                        System.out.println(ex.getMessage());
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }


    //Input: #locus, refineSize
    //Output: double[] ll List for Locus_i
    public double[] getLocusLL(int locusNum, int msSize) throws IOException {
        double[] ll = new double[msSize];
        //String outputFile = _RAxMLdir + "RAxML_result.TEST" + i;
        for(int j = 0; j<msSize; j++){
            String llFile = _RAxMLdir + "RAxML_log.L" + String.valueOf(locusNum) + "T" + String.valueOf(j);
            BufferedReader llReader = new BufferedReader(new FileReader(llFile));
            String[] llString = llReader.readLine().trim().split(" ");
            double tempLL = Double.valueOf(llString[1]);
            ll[j] = tempLL;
            llReader.close();
        }
        return ll;
    }


    //Input: trees[], ll[] ,lociNum
    //Output: double[] ll, trees[]
    public double[] getInitTree(String[] trees, double[] ll, int lociNum) throws IOException {

        for(int j = 0; j<lociNum; j++){
            String tFile = _RAxMLdir + "RAxML_bestTree.T" + String.valueOf(j);
            BufferedReader tReader = new BufferedReader(new FileReader(tFile));
            trees[j] = tReader.readLine().trim();
            tReader.close();
            String llFile = _RAxMLdir + "RAxML_info.T" + String.valueOf(j);
            BufferedReader llReader = new BufferedReader(new FileReader(llFile));
            while(true){
                String temp = llReader.readLine().trim();
                if(temp.startsWith("Final")){
                    String temps[] = temp.split(" ");
                    temp = temps[temps.length-1];
                    ll[j] = Double.valueOf(temp);
                    break;
                }
            }
            llReader.close();
            //System.out.println(trees[j] + "  " + ll[j]);
        }
        return ll;
    }

    //Input: #locus, #bestMSNum
    //Output: RAxML tree
    public Tree getbestGTi(int locusNum, int msTreeNum) throws IOException {
        String gtFile = _RAxMLdir + "RAxML_result.L" + String.valueOf(locusNum) + "T" + String.valueOf(msTreeNum);
        BufferedReader gtReader = new BufferedReader(new FileReader(gtFile));
        String gtString = gtReader.readLine().trim();
        Tree bestGT = Trees.readTree(gtString);
        gtReader.close();
        return bestGT;
    }


    //TODO check
    public static List<UltrametricTree> scaleGTS(List<UltrametricTree> uList, double ratio, boolean undo) {
        Iterator gList = uList.iterator();
        while (gList.hasNext()) {
            UltrametricTree g = (UltrametricTree) gList.next();
            Iterator itNode = g.getTree().getNodes().iterator();
            if (undo == false) {
                while (itNode.hasNext()) {
                    STINode n = (STINode) itNode.next();
                    if (n == (STINode) g.getTree().getRoot())
                        continue;
                    if (n.getParentDistance() != TMutableNode.NO_DISTANCE) {
                        n.setParentDistance(n.getParentDistance() * ratio); // height = height * ratio
                    }
                }
            } else {
                while (itNode.hasNext()) {
                    STINode n = (STINode) itNode.next();
                    if (n == (STINode) g.getTree().getRoot())
                        continue;
                    if (n.getParentDistance() != TMutableNode.NO_DISTANCE) {
                        n.setParentDistance(n.getParentDistance() / ratio); // height = height / ratio
                    }
                }
            }

        }
        return uList;
    }

    public boolean isRooted(Tree t) {
        boolean flag = true;
        TNode root = t.getRoot();
        if(root.getChildCount()>2)
            flag = false;
        Iterator cIt = root.getChildren().iterator();
        TNode c1 = (TNode) cIt.next();
        TNode c2 = (TNode) cIt.next();
        if(c1.getName().equals("O")||c2.getName().equals("O"))
            flag =  true;
        else
            flag = false;
        return flag;
    }

    //scale mtDNA
    public static List<Tree> scaleGTList(boolean multiply, List<Tree> tList, double ratio) {
        List<Tree> newGTS = new ArrayList<Tree>();
        for (int i = 0; i < tList.size(); i++) {
            newGTS.add(scaleGT(tList.get(i), ratio, multiply));
        }
        return tList;
    }

    public List<String> getMSTopos(List<Tree> msGTS, List<String> checkedTopo){
        List<String> topos = new ArrayList<String>();
        for(int i = 0; i<msGTS.size(); i++){
            String temp = getTopo(msGTS.get(i).toString());
            if(checkedTopo.contains(temp))
                continue;
            else
                topos.add(getTopo(msGTS.get(i).toString()));
        }
        return topos;
    }

    public String getTopo(String tree){
        Tree temp = Trees.readTree(tree);
        Iterator node = temp.getNodes().iterator();
        while(node.hasNext()){
            TNode n = (TNode) node.next();
            if(!n.isRoot())
                //TODO 7 26
                n.setParentDistance(Double.NEGATIVE_INFINITY);
        }
        return temp.toString();
    }

    // ./raxmlHPC-PTHREADS -f e -t RAxML_bestTree.T0 -m GTRGAMMA -s 0/dna.phy -n TEST
   // rm *.TEST
//    public static List<Tree> addGTOG(List<Tree> gts, List<Double> ogHeight) {
//        List<Tree> newGTS = new ArrayList<Tree>();
//        for (int i = 0; i < gts.size(); i++) {
//            STITree thisTree = (STITree) gts.get(i);
//            STINode root = thisTree.getRoot();
////            STINode newroot = root.createChild("I" + thisTree.getNodeCount());
////            root.removeChild((TMutableNode) newroot,false);
////            newroot.adoptChild((TMutableNode)root);
////            STINode outgroup = newroot.createChild("O");
////            outgroup.setParentDistance(ogHeight.get(i));
////            root.setParentDistance(Double.NEGATIVE_INFINITY);
////            thisTree.rerootTreeAtNode(newroot);
////            newGTS.add((Tree)thisTree);
//
//            if (thisTree.getNode("O") != null) {
//                newGTS.add((Tree) thisTree);
//                continue;
//            }
//            STINode outgroup = root.createChild("O");
//            outgroup.setParentDistance(ogHeight.get(i));
//            newGTS.add((Tree) thisTree);
//        }
//        return newGTS;
//    }

//    public static List<Tree> removeGTOG(List<Tree> gts) {
//        List<Tree> newGTS = new ArrayList<Tree>();
//        for (int i = 0; i < gts.size(); i++) {
//            MutableTree thisTree = (MutableTree) gts.get(i);
//            TMutableNode root = thisTree.getRoot();
////            TNode outGroup = thisTree.getNode("O");
////            Iterator it = root.getChildren().iterator();
////            TNode newroot = (TNode)it.next();
////            if(newroot.getName().equals("O"))
////                newroot = (TNode)it.next();
////
////            thisTree.rerootTreeAtNode(newroot);
//            root.removeChild(thisTree.getNode("O"), false);
//            newGTS.add(thisTree);
//        }
//        return newGTS;
//    }

    public List<Tree> getOGGTS(List<Tree> gts, double ogHeight){
        List<Tree> ogGTS = new ArrayList<Tree>();
        for(int i = 0;i<gts.size();i++){
            String temp = gts.get(i).toString();
            int length = temp.length();
            Tree ogTree = Trees.readTree("(" + temp.substring(0,length-1) + ",O:" + String.valueOf(ogHeight) + ");");
            ogGTS.add(ogTree);
        }
        return ogGTS;
    }

    //DO NOT change original trees
    public List<Tree> getNoOGGTS(List<Tree> ogGTS){
        List<Tree> gts = new ArrayList<Tree>();
        for(int i = 0;i<ogGTS.size();i++){
            String temp = ogGTS.get(i).toString();
            Tree gt = Trees.readTree(rerootAndRemove(temp, "O"));
            gts.add(gt);
        }
        return gts;
    }

    //Input: rooted binary tree with oG in the first level
    //Output: rooted tree without OG
    public static String rerootAndRemove(String tree, String outgroup) {
        STITree t = (STITree) Trees.readTree(tree);
        if (t == null || t.getNode(outgroup) == null) {
            System.err.println(tree + "  " + outgroup);
            return null;
        }
        t.rerootTreeAtEdge(outgroup);
        t.removeNode(outgroup);
        Trees.removeBinaryNodes(t);
        return t.toNewick();
    }

    //reroot RAxML tree whose root has three children. DO NOT remove OG
    //Input: RAxML tree
    //Output: rooted tree and outgroup is the first-level child
    public static Tree rerootRAxML(Tree t, String outgroup) {
        Trees.autoLabelNodes((MutableTree) t);
        TNode og = t.getNode(outgroup);
        double oHeight = og.getParentDistance();
        String tempP = og.getParent().getName();
        TNode oldRoot = t.getRoot();
        int count = t.getNodeCount();
        if(oldRoot.getChildCount()!=3)
            System.err.println("Special RAxML tree!(2-node root)");
        BniNetwork tempNet = (BniNetwork)Networks.readNetwork(t.toString());
        BniNetNode newRoot = new BniNetNode("I" + String.valueOf(count),0.0);
        BniNetNode oldP = (BniNetNode)tempNet.findNode(tempP);

        oldP.adoptChild(newRoot,oHeight/100);
        newRoot.adoptChild((BniNetNode)tempNet.findNode(outgroup),oHeight);
        oldP.removeChild((BniNetNode)tempNet.findNode(outgroup));

        BniNetNode oldNetRoot = (BniNetNode)tempNet.getRoot();
        BniNetNode tempChild = newRoot;
        BniNetNode tempParent = oldP;
        while(!tempChild.equals(oldNetRoot)){
            Double tempD = tempChild.getParentDistance(tempParent);
            tempParent.removeChild(tempChild);
            tempChild.adoptChild(tempParent,tempD);
            tempChild = tempParent;
            tempParent = (BniNetNode) tempChild.getParents().iterator().next();
        }
        tempNet.resetRoot(newRoot);
        //Tree tempT = tempNet.
        t = Trees.readTree(tempNet.toString());
        return Trees.readTree(tempNet.toString());

    }

    //DO REMOVE OG
    public static Network reRootAST(Tree t, String outgroup) {

        Trees.autoLabelNodes((MutableTree)t);
        TMutableNode oldRoot = (TMutableNode)t.getRoot();
        if(oldRoot.getChildCount()!=2){
            System.err.print("Wrong child number for AST Tree!");
        }
        Iterator it = oldRoot.getChildren().iterator();
        TMutableNode newRoot = oldRoot;
        TMutableNode leafNode = oldRoot;
        while(it.hasNext()){
            TMutableNode temp = (TMutableNode)it.next();
            if(temp.isLeaf())
                leafNode = temp;
            else
                newRoot = temp;
        }
        newRoot.adoptChild(leafNode);
        t.rerootTreeAtNode(newRoot);
        newRoot.removeChild(oldRoot,false);
        Network st = Networks.readNetwork(rerootRAxML(t,outgroup).toNewick());
        NetNode og = st.findNode(outgroup);
        Iterator sibIt = st.getRoot().getChildren().iterator();
        NetNode sib = (NetNode) sibIt.next();
        if(sib.equals(og))
            sib = (NetNode) sibIt.next();
        st.resetRoot(sib);
        return st;
    }



    //UNDO: true -> *halfTheta    ;   false -> /halfTheta
    public static Tree scaleGT(Tree t, double ratio, boolean multiply) {

            Iterator itNode = t.getNodes().iterator();
            if (multiply == false) {
                while (itNode.hasNext()) {
                    STINode n = (STINode) itNode.next();
                    if (n == (STINode) t.getRoot())
                        continue;
                    if (n.getParentDistance() != TMutableNode.NO_DISTANCE) {
                        n.setParentDistance(n.getParentDistance() / ratio); // height = height * ratio
                    }
                }
            } else {
                while (itNode.hasNext()) {
                    STINode n = (STINode) itNode.next();
                    if (n == (STINode) t.getRoot())
                        continue;
                    if (n.getParentDistance() != TMutableNode.NO_DISTANCE) {
                        n.setParentDistance(n.getParentDistance() * ratio); // height = height / ratio
                    }
                }
            }
        return t;
    }

    public List<Tree> getTreeByUTree(List<UltrametricTree> uList) {
        Iterator uts = uList.iterator();
        List<Tree> gts = new ArrayList<Tree>();
        while (uts.hasNext()) {
            STITree gTree = new STITree(((UltrametricTree) uts.next()).getTree());
            gts.add(gTree);
        }
        return gts;
    }

    public List<UltrametricTree> getUTreeByTree(List<Tree> gt) {
        //List<UltrametricTree> ugts = new ArrayList<Tree>();
        Iterator gts = gt.iterator();
        List<UltrametricTree> ultraGTS = new ArrayList<UltrametricTree>();
        while (gts.hasNext()) {
            UltrametricTree uTree = new UltrametricTree((Tree) gts.next());
            ultraGTS.add(uTree);
        }
        return ultraGTS;
    }

    //Input: gt, shortest external length(ST height)
    public Tree getUTree(Tree gt, double base) {
        TNode lowLeaf = gt.getNodes().iterator().next();
        double distance = 0.0;
        Iterator nodeIt = gt.getNodes().iterator();
        while(nodeIt.hasNext()){
            TNode temp = (TNode)nodeIt.next();
            if(temp.isLeaf()){
                double tempDistance = getLeafHeight(temp, gt);
                if(tempDistance>distance){
                    distance = tempDistance;
                    lowLeaf = temp;
                }
            }
            else
                continue;
        }
        double height = distance + base;

        Iterator nodeIt1 = gt.getNodes().iterator();
        while(nodeIt1.hasNext()){
            TNode temp = (TNode)nodeIt1.next();
            if(temp.isLeaf()){
                double addition = height - getLeafHeight(temp, gt);
                temp.setParentDistance(temp.getParentDistance()+addition);

            }
        }
        //ltrametricTree ultraGT = new UltrametricTree(gt);
        return gt;
    }

    public Tree cleanName(Tree t) {
        Iterator node = t.getNodes().iterator();
        while(node.hasNext()){
            STINode stNode = (STINode)node.next();
            if(stNode.isLeaf())
                continue;
            stNode.setName("");
        }
        return t;
    }

    public Network getUTree(Network tempST, double[] height) {
        Tree st = Trees.readTree(tempST.toString());
        TNode lowLeaf = st.getNodes().iterator().next();
        double distance = 0.0;
        Iterator nodeIt = st.getNodes().iterator();
        while(nodeIt.hasNext()){
            TNode temp = (TNode)nodeIt.next();
            if(temp.isLeaf()){
                double tempDistance = getLeafHeight(temp, st);
                if(tempDistance>distance){
                    distance = tempDistance;
                    lowLeaf = temp;
                }
            }
            else
                continue;
        }
        height[2] = distance;

        Iterator nodeIt1 = st.getNodes().iterator();
        while(nodeIt1.hasNext()){
            TNode temp = (TNode)nodeIt1.next();
            if(temp.isLeaf()){
                double addition = height[2] - getLeafHeight(temp, st);
                temp.setParentDistance(temp.getParentDistance()+addition);

            }
        }
        return Networks.readNetwork(st.toString());
    }


    public double getDistance(Tree tree1, Tree tree2) {
        SymmetricDifference symmetricDifference = new SymmetricDifference();
        symmetricDifference.computeDifference(tree1, tree2, false);
        double diff = symmetricDifference.getWeightedAverage();
        return diff;
    }

    public double getRootDistance(Tree tree1, Tree tree2) {
        SymmetricDifference symmetricDifference = new SymmetricDifference();
        symmetricDifference.computeDifference(tree1, tree2, true);
        double diff = symmetricDifference.getWeightedAverage();
        return diff;
    }

    public double sumDoubleList(List<Double> list) {
        double sum = 0.0;
        Iterator it = list.iterator();
        while (it.hasNext()) {
            sum += (double) it.next();
        }
        return sum;
    }

    public double getLeafHeight(NetNode leaf, Network currentST) {
        double currentHeight = 0.0;
        //NetNode leaf = (NetNode) currentST.getLeaves().iterator().next();
        while (leaf != currentST.getRoot()) {
            currentHeight += leaf.getParentDistance((NetNode) leaf.getParents().iterator().next());
            leaf = (NetNode) leaf.getParents().iterator().next();
        }
        return currentHeight;
    }

    //for Tree
    public double getLeafHeight(TNode leaf, Tree currentST) {
        double currentHeight = 0.0;
        //NetNode leaf = (NetNode) currentST.getLeaves().iterator().next();
        while (leaf != currentST.getRoot()) {
            currentHeight += leaf.getParentDistance();
            leaf = (TNode) leaf.getParent();
        }
        return currentHeight;
    }

    public double getNodeHeight(Network currentST, double treeHeight, NetNode<Double> n) {
        //double currentHeight = getNetHeight(currentST);
        double height = 0.0;
        NetNode<Double> leaf = (NetNode<Double>) n.getParents().iterator().next();
        while (leaf != currentST.getRoot()) {
            height += leaf.getParentDistance((NetNode<Double>) leaf.getParents().iterator().next());
            leaf = (NetNode<Double>) leaf.getParents().iterator().next();
        }
        n.setParentDistance((NetNode<Double>) n.getParents().iterator().next(),treeHeight-height);
        return treeHeight-height;
    }

    //
    //
    ////////////////////////////////////
    //For UPRA:

    //Input: MS sample trees
    //Output: distinguished MS topos
    public List<Tree> getDistinguishedTopo(List<Tree> msTrees){
        List<Tree> dTopos = new ArrayList<Tree>();
        List<Integer> sameTopoNum = new ArrayList<Integer>();

        for(int i = 0; i < msTrees.size(); i++){
            if(sameTopoNum.contains(i))
                continue;
            Tree t1 = msTrees.get(i);
            dTopos.add(t1);

            for(int j = i+1 ; j < msTrees.size(); j++){
                if(sameTopoNum.contains(j))
                    continue;
                Tree t2 = msTrees.get(j);
                if(Trees.haveSameRootedTopology(t1,t2)){
                    sameTopoNum.add(j);
                }
            }

        }

        return dTopos;
    }



}
