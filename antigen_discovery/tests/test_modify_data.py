import unittest
import os
import sys
import pandas
import glob
import nepitope
from nepitope import scoring_utils
from nepitope import pep_utils
from nepitope.modify_data import ModifyData


fasta_name = os.path.dirname(os.path.realpath('__file__')) + '/tests/test_fasta.fasta'
swaps_csv_name = os.path.dirname(os.path.realpath('__file__')) + '/tests/summary_results_per_prediction.csv'



class TestModifyData(unittest.TestCase):

    def setUp(self):

        self.test_list = [[0.7815955555555556, 'MDKKYSIGL', 9],
                          [0.7633500000000001, 'DKKYSIGLD', 9],
                          [0.8262666666666667, 'IGTNSVGWA', 9],
                          [0.7322544444444445, 'LEESFLVEE', 9],
                          [0.7556699999999998, 'EESFLVEED', 9],
                          [0.7949977777777778, 'ESFLVEEDK', 9],
                          [0.7851955555555556, 'SFLVEEDKK', 9],
                          [0.7558333333333334, 'FLVEEDKKH', 9],
                          [0.7114655555555555, 'LVEEDKKHE', 9],
                          [0.7133622222222222, 'VEEDKKHER', 9]]

        self.df_cols = ['Align_Col_Number', 'Score', 'Column', 'Peptide']

        self.score_pep_val = [[ 'MDKKYSIGL', len('MDKKYSIGL')]]

        self.swaps_df = pandas.DataFrame([[73, 10, 'nmer', 128, 48, 8, 16, 'RRYTGWGRLS', 653,
        "['RPYTGWGRLS', 'RDYTGWGRLS', 'RFYTGWGRLS', 'RGYTGWGRLS', 'REYTGWGRLS', 'RCYTGWGRLS', 'RIYTGWGRLS', 'RNYTGWGRLS', 'RVYTGWGRLS', 'RAYTGWGRLS', 'RTYTGWGRLS', 'RHYTGWGRLS', 'RSYTGWGRLS', 'RYYTGWGRLS', 'RLYTGWGRLS', 'PRYTGWGRLS', 'RMYTGWGRLS', 'RWYTGWGRLS', 'DRYTGWGRLS', 'RKYTGWGRLS', 'ERYTGWGRLS', 'RQYTGWGRLS', 'NRYTGWGRLS', 'TRYTGWGRLS', 'VRYTGWGRLS', 'RRYTGWGRDS', 'IRYTGWGRLS', 'SRYTGWGRLS', 'RRPTGWGRLS', 'RRDTGWGRLS', 'MRYTGWGRLS', 'ARYTGWGRLS', 'RRETGWGRLS', 'QRYTGWGRLS', 'RRYTGWGRNS', 'LRYTGWGRLS', 'RRYTGWGRES', 'RRYTGDGRLS', 'RRRTGWGRLS', 'GRYTGWGRLS']",
        [653, 654, 655, 656, 657, 658, 659, 660, 661, 662]],
       [52, 9, 'nmer', 102, 45, 22, 11, 'LTFRIPYYV', 443,
        "['LTFRIPYYD', 'LTFRIPYYH', 'LTFRIPYYN', 'LTFRIPYYE', 'LTFRIPYYK', 'LTFRIPYYR', 'LTFRIPYYQ', 'LTFRIPYYY', 'LTFRIPYYW', 'LTFRIPYYP', 'LTFRIPYYG', 'LTFRIPYYF', 'LTFRIPYYS', 'LTFRIPYYC', 'LPFRIPYYV', 'LHFRIPYYV', 'LRFRIPYYV', 'LDFRIPYYV', 'DTFRIPYYV', 'PTFRIPYYV', 'LTRRIPYYV', 'LNFRIPYYV', 'LEFRIPYYV', 'LTKRIPYYV', 'LKFRIPYYV', 'LGFRIPYYV', 'LYFRIPYYV', 'LTFRIPYYT', 'LTFRIPKYV', 'LTERIPYYV', 'LTFRIPRYV', 'ETFRIPYYV', 'LTFRIPGYV', 'LWFRIPYYV', 'LTFRIPYDV', 'LCFRIPYYV']",
        [443, 444, 445, 446, 447, 448, 449, 450, 451]],
       [51, 10, 'nmer', 158, 13, 25, 4, 'LQNEKLYLYY', 805,
        "['LPNEKLYLYY', 'LHNEKLYLYY', 'LRNEKLYLYY', 'LYNEKLYLYY', 'LCNEKLYLYY', 'LKNEKLYLYY', 'LDNEKLYLYY', 'LFNEKLYLYY', 'LNNEKLYLYY', 'DQNEKLYLYY', 'LQNEKLYLYD', 'LQNEKLYLYT', 'LQNEKLYLYE', 'LQNEKLYLYK', 'LQNEKLYLYR', 'LQNEKLYLYN', 'LQNEKLYLYS', 'LQNEKLYLYC', 'LQNEKLYLYP', 'LGNEKLYLYY', 'PQNEKLYLYY', 'LQNEKLYLYQ', 'LQNEKLYLYA', 'LQNEKLYLYV', 'LWNEKLYLYY', 'LQNEKLYLYG', 'LQNEKLYLYH', 'LENEKLYLYY', 'LQNEKLYLYW', 'LQNEKLYLYI', 'LQDEKLYLYY', 'LQEEKLYLYY', 'LANEKLYLYY', 'LQNEKLYLYL', 'LTNEKLYLYY', 'LVNEKLYLYY', 'EQNEKLYLYY', 'LSNEKLYLYY', 'LINEKLYLYY', 'NQNEKLYLYY']",
        [805, 806, 807, 808, 809, 810, 811, 812, 813, 814]],
       [97, 10, 'nmer', 134, 49, 7, 10, 'YLQNGRDMYV', 814,
        "['YLQNGRDMYD', 'YLQNGRDMYH', 'YLQNGRDMYE', 'YLQNGRDMYN', 'YLQNGRDMYK', 'YLQNGRDMYR', 'YLQNGRDMYQ', 'YLQNGRDMYY', 'YLQNGRDMYP', 'YLQNGRDMYW', 'YLQNGRDMYG', 'PLQNGRDMYV', 'DLQNGRDMYV', 'YLQNGRDMYS', 'ELQNGRDMYV', 'YLQNGRDMYC', 'YLQNGRDMYF', 'YLQNGRDMYT', 'YKQNGRDMYV', 'YWQNGRDMYV', 'YHQNGRDMYV', 'YGQNGRDMYV', 'YEQNGRDMYV', 'YDQNGRDMYV', 'YNQNGRDMYV', 'YRQNGRDMYV', 'YPQNGRDMYV', 'YSQNGRDMYV', 'YCQNGRDMYV', 'YYQNGRDMYV', 'YAQNGRDMYV', 'YTQNGRDMYV', 'YFQNGRDMYV', 'QLQNGRDMYV', 'NLQNGRDMYV', 'YLQNGRDKYV', 'HLQNGRDMYV', 'YLQNGRDMDV', 'YVQNGRDMYV', 'CLQNGRDMYV']",
        [814, 815, 816, 817, 818, 819, 820, 821, 822, 823]],
       [80, 10, 'nmer', 175, 11, 7, 7, 'SLLYEYFTVY', 511,
        "['SLLYEYFTVD', 'SLLYEYFTVR', 'SLLYEYFTVK', 'SLLYEYFTVN', 'SLLYEYFTVE', 'SLLYEYFTVT', 'SLLYEYFTVS', 'SLLYEYFTVP', 'SLLYEYFTVC', 'SLLYEYFTVQ', 'SLLYEYFTVA', 'SLLYEYFTVV', 'SLLYEYFTVG', 'SLLYEYFTVH', 'SLLYEYFTVW', 'SLLYEYFTVI', 'SLLYEYFTWY', 'SLLYEYFTVL', 'SLLYEYFTKY', 'SLLYEYFTDY', 'SLLYEYFTLY', 'SLLYEYKTVY', 'SLLYEYFTRY', 'SLLYEYETVY', 'SLLYEYFTMY', 'SLLYEYPTVY', 'SLLYEYFTQY', 'SLLYEYQTVY', 'SLLYEYFTVM', 'SLLYEYFTIY', 'SLLYEYRTVY', 'SLLLEYFTVY', 'SLLEEYFTVY', 'SLLYEYFWVY', 'SLLYEYFTFY', 'SLLYEYFTGY', 'SLLYEYATVY', 'SLLYEYFTNY', 'SLLYEYHTVY', 'SLLTEYFTVY']",
        [511, 512, 513, 514, 515, 516, 517, 518, 519, 520]]], columns = ['Unnamed: 0', 'nmer', 'allele', 'num high affinity peps',
       'num med affinity peps', 'num low affinity peps',
       'num no affinity peps', 'original peptide', 'original pos',
       'top scoring peptides', 'Range'])


        self.test_fasta = fasta_name
        self.csv = swaps_csv_name
        self.hashed_fasta = [['S__pyogenes_Cas9', 'MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGALLFDSGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFHRLEESFLVEEDKKHERHPIFGNIVDEVAYHEKYPTIYHLRKKLVDSTDKADLRLIYLALAHMIKFRGHFLIEGDLNPDNSDVDKLFIQLVQTYNQLFEENPINASGVDAKAILSARLSKSRRLENLIAQLPGEKKNGLFGNLIALSLGLTPNFKSNFDLAEDAKLQLSKDTYDDDLDNLLAQIGDQYADLFLAAKNLSDAILLSDILRVNTEITKAPLSASMIKRYDEHHQDLTLLKALVRQQLPEKYKEIFFDQSKNGYAGYIDGGASQEEFYKFIKPILEKMDGTEELLVKLNREDLLRKQRTFDNGSIPHQIHLGELHAILRRQEDFYPFLKDNREKIEKILTFRIPYYVGPLARGNSRFAWMTRKSEETITPWNFEEVVDKGASAQSFIERMTNFDKNLPNEKVLPKHSLLYEYFTVYNELTKVKYVTEGMRKPAFLSGEQKKAIVDLLFKTNRKVTVKQLKEDYFKKIECFDSVEISGVEDRFNASLGTYHDLLKIIKDKDFLDNEENEDILEDIVLTLTLFEDREMIEERLKTYAHLFDDKVMKQLKRRRYTGWGRLSRKLINGIRDKQSGKTILDFLKSDGFANRNFMQLIHDDSLTFKEDIQKAQVSGQGDSLHEHIANLAGSPAIKKGILQTVKVVDELVKVMGRHKPENIVIEMARENQTTQKGQKNSRERMKRIEEGIKELGSQILKEHPVENTQLQNEKLYLYYLQNGRDMYVDQELDINRLSDYDVDHIVPQSFLKDDSIDNKVLTRSDKNRGKSDNVPSEEVVKKMKNYWRQLLNAKLITQRKFDNLTKAERGGLSELDKAGFIKRQLVETRQITKHVAQILDSRMNTKYDENDKLIREVKVITLKSKLVSDFRKDFQFYKVREINNYHHAHDAYLNAVVGTALIKKYPKLESEFVYGDYKVYDVRKMIAKSEQEIGKATAKYFFYSNIMNFFKTEITLANGEIRKRPLIETNGETGEIVWDKGRDFATVRKVLSMPQVNIVKKTEVQTGGFSKESILPKRNSDKLIARKKDWDPKKYGGFDSPTVAYSVLVVAKVEKGKSKKLKSVKELLGITIMERSSFEKNPIDFLEAKGYKEVKKDLIIKLPKYSLFELENGRKRMLASAGELQKGNELALPSKYVNFLYLASHYEKLKGSPEDNEQKQLFVEQHKHYLDEIIEQISEFSKRVILADANLDKVLSAYNKHRDKPIREQAENIIHLFTLTNLGAPAAFKYFDTTIDRKRYTSTKEVLDATLIHQSITGLYETRIDLSQLGGD'],
                             ['Staphylococcus_aureus_Cas9', 'MKRNYILGLDIGITSVGYGIIDYETRDVIDAGVRLFKEANVENNEGRRSKRGARRLKRRRRHRIQRVKKLLFDYNLLTDHSELSGINPYEARVKGLSQKLSEEEFSAALLHLAKRRGVHNVNEVEEDTGNELSTKEQISRNSKALEEKYVAELQLERLKKDGEVRGSINRFKTSDYVKEAKQLLKVQKAYHQLDQSFIDTYIDLLETRRTYYEGPGEGSPFGWKDIKEWYEMLMGHCTYFPEELRSVKYAYNADLYNALNDLNNLVITRDENEKLEYYEKFQIIENVFKQKKKPTLKQIAKEILVNEEDIKGYRVTSTGKPEFTNLKVYHDIKDITARKEIIENAELLDQIAKILTIYQSSEDIQEELTNLNSELTQEEIEQISNLKGYTGTHNLSLKAINLILDELWHTNDNQIAIFNRLKLVPKKVDLSQQKEIPTTLVDDFILSPVVKRSFIQSIKVINAIIKKYGLPNDIIIELAREKNSKDAQKMINEMQKRNRQTNERIEEIIRTTGKENAKYLIEKIKLHDMQEGKCLYSLEAIPLEDLLNNPFNYEVDHIIPRSVSFDNSFNNKVLVKQEENSKKGNRTPFQYLSSSDSKISYETFKKHILNLAKGKGRISKTKKEYLLEERDINRFSVQKDFINRNLVDTRYATRGLMNLLRSYFRVNNLDVKVKSINGGFTSFLRRKWKFKKERNKGYKHHAEDALIIANADFIFKEWKKLDKAKKVMENQMFEEKQAESMPEIETEQEYKEIFITPHQIKHIKDFKDYKYSHRVDKKPNRELINDTLYSTRKDDKGNTLIVNNLNGLYDKDNDKLKKLINKSPEKLLMYHHDPQTYQKLKLIMEQYGDEKNPLYKYYEETGNYLTKYSKKDNGPVIKKIKYYGNKLNAHLDITDDYPNSRNKVVKLSLKPYRFDVYLDNGVYKFVTVKNLDVIKKENYYEVNSKCYEEAKKLKKISNQAEFIASFYNNDLIKINGELYRVIGVNNDLLNRIEVNMIDITYREYLENMNDKRPPRIIKTIASKTQSIKKYSTDILGNLYEVKSKKHPQIIKKG'],
                             ['S_CRISPR_1_thermophilus_Cas9', 'MSDLVLGLDIGIGSVGVGILNKVTGEIIHKNSRIFPAAQAENNLVRRTNRQGRRLARRKKHRRVRLNRLFEESGLITDFTKISINLNPYQLRVKGLTDELSNEELFIALKNMVKHRGISYLDDASDDGNSSVGDYAQIVKENSKQLETKTPGQIQLERYQTYGQLRGDFTVEKDGKKHRLINVFPTSAYRSEALRILQTQQEFNPQITDEFINRYLEILTGKRKYYHGPGNEKSRTDYGRYRTSGETLDNIFGILIGKCTFYPDEFRAAKASYTAQEFNLLNDLNNLTVPTETKKLSKEQKNQIINYVKNEKAMGPAKLFKYIAKLLSCDVADIKGYRIDKSGKAEIHTFEAYRKMKTLETLDIEQMDRETLDKLAYVLTLNTEREGIQEALEHEFADGSFSQKQVDELVQFRKANSSIFGKGWHNFSVKLMMELIPELYETSEEQMTILTRLGKQKTTSSSNKTKYIDEKLLTEEIYNPVVAKSVRQAIKIVNAAIKEYGDFDNIVIEMARETNEDDEKKAIQKIQKANKDEKDAAMLKAANQYNGKAELPHSVFHGHKQLATKIRLWHQQGERCLYTGKTISIHDLINNSNQFEVDHILPLSITFDDSLANKVLVYATANQEKGQRTPYQALDSMDDAWSFRELKAFVRESKTLSNKKKEYLLTEEDISKFDVRKKFIERNLVDTRYASRVVLNALQEHFRAHKIDTKVSVVRGQFTSQLRRHWGIEKTRDTYHHHAVDALIIAASSQLNLWKKQKNTLVSYSEDQLLDIETGELISDDEYKESVFKAPYQHFVDTLKSKEFEDSILFSYQVDSKFNRKISDATIYATRQAKVGKDKADETYVLGKIKDIYTQDGYDAFMKIYKKDKSKFLMYRHDPQTFEKVIEPILENYPNKQINDKGKEVPCNPFLKYKEEHGYIRKYSKKGNGPEIKSLKYYDSKLGNHIDITPKDSNNKVVLQSVSPWRADVYFNKTTGKYEILGLKYADLQFDKGTGTYKISQEKYNDIKKKEGVDSDSEFKFTLYKNDLLLVKDTETKEQQLFRFLSRTMPKQKHYVELKPYDKQKFEGGEALIKVLGNVANSGQCKKGLGKSNISIYKVRTDVLGNQHIIKNEGDKPKLDF'],
                             ['N__meningitidis_Cas9', 'MAAFKPNPINYILGLDIGIASVGWAMVEIDEDENPICLIDLGVRVFERAEVPKTGDSLAMARRLARSVRRLTRRRAHRLLRARRLLKREGVLQAADFDENGLIKSLPNTPWQLRAAALDRKLTPLEWSAVLLHLIKHRGYLSQRKNEGETADKELGALLKGVADNAHALQTGDFRTPAELALNKFEKESGHIRNQRGDYSHTFSRKDLQAELILLFEKQKEFGNPHVSGGLKEGIETLLMTQRPALSGDAVQKMLGHCTFEPAEPKAAKNTYTAERFIWLTKLNNLRILEQGSERPLTDTERATLMDEPYRKSKLTYAQARKLLGLEDTAFFKGLRYGKDNAEASTLMEMKAYHAISRALEKEGLKDKKSPLNLSPELQDEIGTAFSLFKTDEDITGRLKDRIQPEILEALLKHISFDKFVQISLKALRRIVPLMEQGKRYDEACAEIYGDHYGKKNTEEKIYLPPIPADEIRNPVVLRALSQARKVINGVVRRYGSPARIHIETAREVGKSFKDRKEIEKRQEENRKDREKAAAKFREYFPNFVGEPKSKDILKLRLYEQQHGKCLYSGKEINLGRLNEKGYVEIDHALPFSRTWDDSFNNKVLVLGSENQNKGNQTPYEYFNGKDNSREWQEFKARVETSRFPRSKKQRILLQKFDEDGFKERNLNDTRYVNRFLCQFVADRMRLTGKGKKRVFASNGQITNLLRGFWGLRKVRAENDRHHALDAVVVACSTVAMQQKITRFVRYKEMNAFDGKTIDKETGEVLHQKTHFPQPWEFFAQEVMIRVFGKPDGKPEFEEADTPEKLRTLLAEKLSSRPEAVHEYVTPLFVSRAPNRKMSGQGHMETVKSAKRLDEGVSVLRVPLTQLKLKDLEKMVNREREPKLYEALKARLEAHKDDPAKAFAEPFYKYDKAGNRTQQVKAVRVEQVQKTGVWVRNHNGIADNATMVRVDVFEKGDKYYLVPIYSWQVAKGILPDRAVVQGKDEEDWQLIDDSFNFKFSLHPNDLVEVITKKARMFGYFASCHRGTGNINIRIHDLDHKIGKNGILEGIGVKTALSFQKYQIDELGKEIRPCRLKKRPPVR']]


    def test_scoring_utils(self):
        self.assertEqual(self, 1, 2)

    def test_hash_dic(self):
        swap_dic = ModifyData._hash_swaps()
        self.assertNotEqual(list(swap_dic.values()), list(swap_dic.keys()))
        self.assertEqual([len(i) for i in swap_dic.values()], [len(i) for i in swap_dic.keys()])

    def test_hash_fasta(self):
        swap_dic = ModifyData._hash_fasta()
        self.assertEqual()






