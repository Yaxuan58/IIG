package edu.rice.cs.bioinfo.programs.phylonet.algos.iterHeuristic;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by doriswang on 1/23/18.
 */
public class Treefix {

        private  String OutDir;
        private  String InputDir;
        private  String RAxML_dir;

        private  String rTree_DIR;
        private int TREEINDEX;




        Treefix(){
//
//        private Network<NetNodeInfo> INIT_ST; //cUnit
//        private Network<NetNodeInfo> currentST;
//        private static String trueST;
//        private Tree Constraint_ST;// GLASS tree from UPGMA. For external branch
//
//        private List<String> trueGTS;
//        private List<Tree> trueTopos;
//        private List<Tree> finalGTS;
//        private static List<Alignment> trueSeq;
//
//
//        private String GLOBALBESTST;
//        private List<String> GLOBALBESTGTS;
    }
//    public void runUpdateShell(int msSize, int lociNum) throws IOException {
//
//        //String seqPath = _RAxMLdir;
//        for (int i = 0; i < lociNum; i++) {
//            String seqPath = _RAxMLdir  + i + "/";
//            for(int j = 0; j<msSize; j++){
//                //updateL0T1.sh
//                String cmdFile = _RAxMLdir  + i + "/updateL" + String.valueOf(i) + "T" + String.valueOf(j) + ".sh";
//                ProcessBuilder pb = new ProcessBuilder("/bin/bash", cmdFile);
//                pb.redirectErrorStream(true);
//                try {
//                    Process proc = pb.start();
//                    try {
//                        proc.waitFor();
//                    } catch (InterruptedException ex) {
//                        System.out.println(ex.getMessage());
//                    }
//                } catch (Exception e) {
//                    e.printStackTrace();
//                }
//            }
//        }
//    }
//    //Input: trees[], ll[] ,lociNum
//    //Output: double[] ll, trees[]
//    public double[] getInitTree(String[] trees, double[] ll, int lociNum) throws IOException {
//
//        for(int j = 0; j<lociNum; j++){
//            String tFile = _RAxMLdir  + "RAxML_bestTree.T" + String.valueOf(j);
//            BufferedReader tReader = new BufferedReader(new FileReader(tFile));
//            trees[j] = tReader.readLine().trim();
//            tReader.close();
//            String llFile = _RAxMLdir + "RAxML_info.T" + String.valueOf(j);
//            BufferedReader llReader = new BufferedReader(new FileReader(llFile));
//            while(true){
//                String temp = llReader.readLine().trim();
//                if(temp.startsWith("Final")){
//                    String temps[] = temp.split(" ");
//                    temp = temps[temps.length-1];
//                    ll[j] = Double.valueOf(temp);
//                    break;
//                }
//            }
//            llReader.close();
//            //System.out.println(trees[j] + "  " + ll[j]);
//        }
//        return ll;
//    }
//    //./raxmlHPC -M -m GTRGAMMA -p 12345 -q simpleDNApartition.txt -s dna.phy -n T22
//    public List<Tree> initByRAxML(List<Alignment> aln, int seqLen, int taxaNum, int lociNum) throws IOException, ParseException, InterruptedException {
//        //writeSeq + partition     rm + raxml    read in
//        List<Tree> gts = new ArrayList<Tree>();
//        for (int i = 0; i < lociNum; i++) {
//            System.out.println("RAXML No" + i + "Start...");
//            isExitsPath(_RAxMLdir + i);
//            String rm = "rm *.T" + i + '\n';
//            //./raxmlHPC-PTHREADS -M -m GTRGAMMA -p 12345 -# 10 -s 0/dna.phy -n T0
//            String raxml = _RAxMLdir + "raxmlHPC-PTHREADS -M -m GTRGAMMA -p 12345 -# 30 " + "-s " + i + "/dna.phy -n T" + i;
//            String cmdFile = _RAxMLdir + i + "/tempCMD.sh";
//            BufferedWriter cmd = new BufferedWriter(new FileWriter(cmdFile));
//            String tempR = _RAxMLdir ;
//            cmd.write("cd " + tempR.substring(0, tempR.length() - 1) + '\n');
//            cmd.write(rm);
//            cmd.flush();
//            cmd.write(raxml);
//            cmd.flush();
//            cmd.close();
//
//            ProcessBuilder pb = new ProcessBuilder("/bin/bash", cmdFile);
//            pb.redirectErrorStream(true);
//            try {
//                Process proc = pb.start();
//                try {
//                    proc.waitFor();
//                } catch (InterruptedException ex) {
//                    System.out.println(ex.getMessage());
//                }
//            } catch (Exception e) {
//                e.printStackTrace();
//            }
//        }
//        for (int i = 0; i < lociNum; i++) {
//            String outputFile = _RAxMLdir + "RAxML_bestTree.T" + i;
//            BufferedReader gtReader = new BufferedReader(new FileReader(outputFile));
//            String gtString = gtReader.readLine().trim();
//            gts.add(Trees.readTree(gtString));
//            gtReader.close();
//        }
//        return gts;
//    }
//    public String[] getRAxMLStartT(int lociNum) throws IOException, ParseException, InterruptedException {
//        //List<Tree> gts = new ArrayList<Tree>();
//        String[] gts = new String[lociNum];
//        for (int i = 0; i < lociNum; i++) {
//            String tree = _RAxMLdir + "/RAxML_Init/RAxML_bestTree.T" + i;
//            String llFile = _RAxMLdir + "/RAxML_Init/RAxML_info.T" + i;
//            BufferedReader gtReader = new BufferedReader(new FileReader(tree));
//            String gtString = gtReader.readLine().trim();
//            gtReader.close();
//            BufferedReader llReader = new BufferedReader(new FileReader(llFile));
//            double ll = 0.0;
//            String llString = "";
//            while(llString!=null){
//                llString = llReader.readLine();
//                if(llString.length()<20||llString.isEmpty())
//                    continue;
//                if(llString.substring(0,5).equals("Final")){
//                    String[] result = llString.split(" ");
//                    ll = Double.valueOf(result[6]);
//                    break;
//                }
//            }
//            gts[i] = gtString + " " + ll;
//            llReader.close();
//        }
//        return gts;
//    }
//    public void getUpdateSh(int lociNum, int msSize) throws IOException, ParseException, InterruptedException {
//        //writeSeq + partition     rm + raxml    read in
//        List<Tree> gts = new ArrayList<Tree>();
//        for (int i = 0; i < lociNum; i++) {
//            String seqPath;
//            for(int j = 0; j<msSize; j++){
//                seqPath = _RAxMLdir  + i + "/";
//                String tName = "L" + String.valueOf(i) + "T" + String.valueOf(j); // L_i T_j
//                BufferedWriter phy = new BufferedWriter(new FileWriter(seqPath + "updateL" + String.valueOf(i) + "T" + String.valueOf(j) + ".sh"));
//                System.out.println(seqPath + "updateL" + String.valueOf(i) + "T" + String.valueOf(j) + ".sh");
//                phy.write("cd " + _RAxMLdir + '\n');
//                phy.flush();
//                phy.write("rm *." + tName + '\n');
//                phy.flush();
//                phy.write(_RAxMLdir    + "raxmlHPC-PTHREADS  -f e -t topo.T" + String.valueOf(j) + " -m GTRGAMMA -s " + String.valueOf(i) + "/dna.phy -n " + tName + '\n');
//                phy.flush();
//                phy.close();
//            }
//
//        }
//        //return gts;
//    }
//    //write alignment with 16 taxa for RAxML
//    public List<Tree> get16Aln(List<Alignment> aln, int seqLen, int taxaNum, int lociNum) throws IOException, ParseException, InterruptedException {
//        //writeSeq + partition     rm + raxml    read in
//        List<Tree> gts = new ArrayList<Tree>();
//        for (int i = 0; i < lociNum; i++) {
//            isExitsPath(_RAxMLdir +i);
//            BufferedWriter phy = new BufferedWriter(new FileWriter(_RAxMLdir  + i + "/dnaNoOG.phy"));
//            phy.write(taxaNum + " " + seqLen + '\n');
//            phy.flush();
//            Alignment a = aln.get(i);
//            List<String> names = a.getTaxaNames();
//            Map<String, String> thisAln = a.getAlignment();
//            for (int j = 0; j < taxaNum; j++) {
//                String name = names.get(j);
//                if (name == "O")
//                    continue;
//                String seq = thisAln.get(name);
//                phy.write(name + '\t' + seq + '\n');
//            }
//            phy.flush();
//            phy.close();
//        }
//        return gts;
//    }

}
