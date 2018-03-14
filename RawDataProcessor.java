package edu.rice.cs.bioinfo.programs.phylonet.algos.iterHeuristic;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.util.HashMap;
import java.util.Iterator;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;

/**
 * Created by doriswang on 9/1/17.
 * process Rattle snack data
 * //RawDataProcessor r = new RawDataProcessor();
   //r.compare();
 */
public class RawDataProcessor {


    public HashMap<Integer, String> getTaxaNameMap() throws IOException{
        HashMap<Integer, String> nameMap = new HashMap<Integer, String>();
        BufferedReader mapfile = new BufferedReader(new FileReader("/Users/doriswang/PhyloNet/Data/syr/" + "mapping_Taxa.txt"));
        for(int i = 0; i<26; i++){
            String line = mapfile.readLine();
            String[] lines = line.split(":");
            nameMap.put(Integer.valueOf(lines[0]),lines[1]);
        }
        return nameMap;
    }

    public String getSTString() throws IOException{
        //String st = "";
        String raw = "(Agkistrodon_contortrix_1:0.03798,Agkistrodon_piscivorus_1:0.03798,(((((Sistrurus_catenatus_catenatus_WI_1:0.004797,Sistrurus_catenatus_catenatus_IL2_1:0.003204,Sistrurus_catenatus_catenatus_IL1_1:0.000806,Sistrurus_catenatus_catenatus_MI_1:0.000772,Sistrurus_catenatus_catenatus_ON1_1:0.000810,Sistrurus_catenatus_catenatus_PA_1:0.000787,Sistrurus_catenatus_catenatus_NY_1:0.000786,Sistrurus_catenatus_catenatus_OH_1:0.000798)0.56:0.001607,Sistrurus_catenatus_catenatus_ON2_1:0.002159)1.00:0.037802,(((Sistrurus_catenatus_edwardsii_NM1_1:0.000817,Sistrurus_catenatus_edwardsii_NM2_1:0.000808,Sistrurus_catenatus_edwardsii_AZ_1:0.000809)1.00:0.003299,(Sistrurus_catenatus_edwardsii_CO_1:0.000806,Sistrurus_catenatus_tergeminus_KS3_1:0.000804,Sistrurus_catenatus_tergeminus_KS2_1:0.000819,Sistrurus_catenatus_tergeminus_KS1_1:0.000806)0.98:0.002343)0.90:0.002431,(Sistrurus_catenatus_tergeminus_MO1_1:0.000788,Sistrurus_catenatus_tergeminus_MO2_1:0.000807)0.88:0.001933)1.00:0.029201)1.00:0.032477,((Sistrurus_miliarius_miliarius_NC_1:0.000816,(Sistrurus_miliarius_barbouri_FL2_1:0.000801,Sistrurus_miliarius_barbouri_FL3_1:0.000811,Sistrurus_miliarius_barbouri_FL1_1:0.000801)1.00:0.004790,Sistrurus_miliarius_streckeri_OK1_1:0.000806):0.001593,Sistrurus_miliarius_streckeri_OK2_1:0.001260)1.00:0.049365)1.00:0.076642)1.00:0.031080);";
        HashMap name = getTaxaNameMap();
        for(int i = 0; i < name.size(); i++){
            String temp = (String)name.get(i);
            int index = raw.indexOf(temp);
            String newRaw = raw.substring(0,index) + String.valueOf(i) + raw.substring(index+temp.length());
            raw = newRaw;
        }
        return raw;
    }

    public void compare() throws IOException{
        String syrPath = "/Users/doriswang/PhyloNet/Data/syr/syr011.txt";
        BufferedReader stReader = new BufferedReader(new FileReader(syrPath));
        for (int ln = 0; ln < 6 ; ln++) {
            String tree = stReader.readLine().trim();
        }
        BufferedReader r1 = new BufferedReader(new FileReader("/Users/doriswang/PhyloNet/Data/syr/1/" + "seq.txt"));
        BufferedReader r2 = new BufferedReader(new FileReader("/Users/doriswang/PhyloNet/Data/syr/2/" + "seq.txt"));
        //BufferedWriter mapping = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/Data/syr/" + "mapping_Taxa.txt"));

        for(int ln = 0; ln < 26 ; ln++) {
            String[] s = r1.readLine().trim().split(":");
            String s1 = s[1];
            s = r2.readLine().trim().split(":");
            String s2 = s[1];
            if(s1 == s2)
                System.out.println(String.valueOf(ln) + "Same");
            else {
                System.out.println("#" + String.valueOf(ln) + " Different: ");
                int count = 0;
                for(int i = 0; i<s1.length(); i++){
                    if(s1.charAt(i)!=s2.charAt(i)) {
                        count++;
                        System.out.print(String.valueOf(i) + " ");
                    }
                }
                System.out.println(" ");
                System.out.println(count + "  differences!");
            }
        }
        r1.close();
        r2.close();
    }

    public void syrRaw() throws IOException{
        String syrPath = "/Users/doriswang/PhyloNet/Data/syr/syr011.txt";
        BufferedReader stReader = new BufferedReader(new FileReader(syrPath));
        for (int ln = 0; ln < 6 ; ln++) {
            String tree = stReader.readLine().trim();
        }
        BufferedWriter writer1 = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/Data/syr/1/" + "syeq1.txt"));
        BufferedWriter writer2 = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/Data/syr/2/" + "seq2.txt"));
        BufferedWriter mapping = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/Data/syr/" + "mapping_Taxa.txt"));

        for(int ln = 0; ln < 52 ; ln++) {
            String seq = stReader.readLine().trim();
            String[] temp = seq.split("\t");
            if(temp[0].substring(temp[0].length()-2).equals("_1")){
                writer1.write(String.valueOf(ln/2) + ":" + temp[1]+"\n");
                mapping.write(String.valueOf(ln/2) + ":" + temp[0] + "\n");
            }
            if(temp[0].substring(temp[0].length()-2).equals("_2")){
                writer2.write(String.valueOf(ln/2) + ":" + temp[1]+"\n");
            }
            //trueGTS.add(tree);
        }

        writer1.close();
        writer2.close();
        mapping.close();
        stReader.close();
    }

    //Input: SYR data
    //Output:seq for RaxML without taxa that misswhole locus data
    public void syrGetSeq() throws IOException{
        String syrPath = "/Users/doriswang/PhyloNet/Data/SYRData/";

        for (int ln = 0; ln < 19 ; ln++) {
            BufferedReader stReader = new BufferedReader(new FileReader(syrPath + ln+"/seq1.phy"));
            BufferedWriter stWriter = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/tools/standard-RAxML-master/" + ln +  "/dna.phy"));

            String info = stReader.readLine().trim();
            info = "21" + info.substring(2);
            int j = 0;
            while(j<26) {
                String temp = info.split("\\t")[0];
                if(temp.equals("0")||temp.equals("3")||temp.equals("15")||temp.equals("21")||temp.equals("11")){
                    j++;
                    info = stReader.readLine().trim();
                    continue;
                }

                stWriter.write(info + "\n");
                stWriter.flush();
                info = stReader.readLine().trim();
                j++;
            }
            stWriter.write(info + "\n");
            stWriter.flush();
            stReader.close();
            stWriter.close();
        }
        for (int ln = 0; ln < 19 ; ln++) {
            BufferedReader stReader = new BufferedReader(new FileReader(syrPath + ln+"/seq2.phy"));
            BufferedWriter stWriter = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/tools/standard-RAxML-master/" + (ln+19) +  "/dna.phy"));

            String info = stReader.readLine().trim();
            info = "21" + info.substring(2);
            int j = 0;
            while(j<26) {
                String temp = info.split("\\t")[0];
                if(temp.equals("0")||temp.equals("3")||temp.equals("15")||temp.equals("21")||temp.equals("11")){
                    j++;
                    info = stReader.readLine().trim();
                    continue;
                }

                stWriter.write(info + "\n");
                stWriter.flush();
                info = stReader.readLine().trim();
                j++;
            }
            stWriter.write(info + "\n");
            stWriter.flush();
            stReader.close();
            stWriter.close();
        }
    }

    public String addTaxonName(String tempStr) throws IOException {
        //String newTree;
        String syrPath = "/Users/doriswang/PhyloNet/Data/SYRData/";
        String mapPath = "/Users/doriswang/PhyloNet/Data/syr/20OGMapping.txt";
        BufferedReader mapReader = new BufferedReader(new FileReader(mapPath));
        HashMap<String, String> tMap = new HashMap<String, String>();
        for (int i = 0; i < 21; i++) {
            //String line =
            String[] lines = mapReader.readLine().trim().split(":");
            String i1 = lines[1];
            String i2 = lines[2];
            tMap.put(i1, i2);
        }
        MutableTree tempTree = (MutableTree) Trees.readTree(tempStr);
        Iterator<TMutableNode> it = (Iterator<TMutableNode>) tempTree.getNodes().iterator();
        while(it.hasNext()){
            TMutableNode n = it.next();
            if(n.isLeaf()) {
                String name = tMap.get(n.getName());
                n.setName(name);
            }
            else
                continue;
        }
        return tempTree.toString();
    }

    public void doubleMap() throws IOException{
        String syrPath = "/Users/doriswang/PhyloNet/Data/SYRData/";
        String mapPath = "/Users/doriswang/PhyloNet/Data/syr/20OGMapping.txt";
        BufferedReader mapReader = new  BufferedReader(new FileReader(mapPath));
        HashMap<String,String> tMap = new HashMap<String,String>();
        for(int i = 0;i< 21;i++){
            //String line =
            String[] lines = mapReader.readLine().trim().split(":");
            String i1 = lines[0];
            String i2 = lines[1];
            tMap.put(i1,i2);
        }
        for (int ln = 0; ln < 19 ; ln++) {
            BufferedReader stReader = new BufferedReader(new FileReader(syrPath + ln+"/seq1.phy"));
            BufferedWriter stWriter = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/tools/standard-RAxML-master/" + ln +  "/dna.phy"));

            String info = stReader.readLine().trim();
            info = "21" + info.substring(2);
            stWriter.write(info + "\n");
            stWriter.flush();
            int j = 0;
            while(j<26) {
                String[] temp = info.split("\\t");
                if (tMap.containsKey(temp[0])) {
                    stWriter.write(tMap.get(temp[0]) + " " + temp[1] + "\n");
                    stWriter.flush();
                    info = stReader.readLine().trim();
                    j++;
                } else {
                    j++;
                    info = stReader.readLine().trim();
                    continue;
                }
            }
            String[] temp = info.split("\\t");
            if (tMap.containsKey(temp[0])) {
                stWriter.write(tMap.get(temp[0]) + " " + temp[1] + "\n");
                stWriter.flush();

            }
        }
        for (int ln = 0; ln < 19 ; ln++) {
            BufferedReader stReader = new BufferedReader(new FileReader(syrPath + ln+"/seq2.phy"));
            BufferedWriter stWriter = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/tools/standard-RAxML-master/" + (ln+19) +  "/dna.phy"));

            String info = stReader.readLine().trim();
            info = "21" + info.substring(2);
            stWriter.write(info + "\n");
            stWriter.flush();
            int j = 0;
            while(j<26) {
                String[] temp = info.split("\\t");
                if (tMap.containsKey(temp[0])) {
                    stWriter.write(tMap.get(temp[0]) + " " + temp[1] + "\n");
                    stWriter.flush();
                    info = stReader.readLine().trim();
                    j++;
                } else {
                    j++;
                    info = stReader.readLine().trim();
                    continue;
                }


            }

            String[] temp = info.split("\\t");
            if (tMap.containsKey(temp[0])) {
                stWriter.write(tMap.get(temp[0]) + " " + temp[1] + "\n");
                stWriter.flush();

            }
            stReader.close();
            stWriter.close();
        }
    }

    //get seq length and name
    public void getSeqInfo() throws IOException{
        String syrPath = "/Users/doriswa 2ng/PhyloNet/Data/syr/syr011.txt";
        BufferedReader stReader = new BufferedReader(new FileReader(syrPath));
        for (int ln = 0; ln < 6 ; ln++) {
            String tree = stReader.readLine().trim();
        }
        BufferedWriter writer1 = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/Data/syr/1/" + "syeq1.txt"));
        BufferedWriter writer2 = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/Data/syr/2/" + "seq2.txt"));
        BufferedWriter mapping = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/Data/syr/" + "mapping_Taxa.txt"));

        for(int ln = 0; ln < 52 ; ln++) {
            String seq = stReader.readLine().trim();
            String[] temp = seq.split("\t");
            if(temp[0].substring(temp[0].length()-2).equals("_1")){
                writer1.write(String.valueOf(ln/2) + ":" + temp[1]+"\n");
                mapping.write(String.valueOf(ln/2) + ":" + temp[0] + "\n");
            }
            if(temp[0].substring(temp[0].length()-2).equals("_2")){
                writer2.write(String.valueOf(ln/2) + ":" + temp[1]+"\n");
            }
            //trueGTS.add(tree);
        }

        writer1.close();
        writer2.close();
        mapping.close();
        stReader.close();
    }
    //TODO 9-4 split different seq   ;  put into different directory which named by locus name

    public void syrSplit() throws  IOException{
        String syrPath = "/Users/doriswang/PhyloNet/Data/syr/";

        BufferedReader r1 = new BufferedReader(new FileReader(syrPath+ "1/seq1.txt"));
        BufferedReader r2 = new BufferedReader(new FileReader(syrPath + "2/seq2.txt"));
        BufferedWriter w1 = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/Data/syr/ATP/1/" + "seq.txt"));
        BufferedWriter w2 = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/Data/syr/ATP/2/" + "seq.txt"));
        RawDataProcessor rdp = new RawDataProcessor();
        for(int ln = 0; ln < 26 ; ln++) {
            String seq = r1.readLine().trim();
            String[] temp = seq.split(":");
            if(temp[1].length()<2)
                temp[1] = temp[2];
            w1.write(temp[0]+  ":" + temp[1].substring(0,665)+"\n");
            seq = r2.readLine().trim();
            temp = seq.split(":");
            if(temp[1].length()<2)
                temp[1] = temp[2];
            w2.write(temp[0]+  ":" + temp[1].substring(0,665)+"\n");
        }
        w1.close();
        w2.close();
        r1.close();
        r2.close();
    }
}
