%% Compare periCSD ROI-level fluor signals
% Identify ROI within an axon that are suited for measuring conduction velocity. Order from left to right
% Construct axon fibers 
fiber = cell(1,Nexpt);
%DL102 180426
fiber{1}(1).ROI = [48,46,42]; 
fiber{1}(2).ROI = [17,31,36,45,41]; 
fiber{1}(3).ROI = [16,25,23,26,14,32]; 
fiber{1}(4).ROI = [10,13,15]; 
fiber{1}(5).ROI = [7,8,6,11]; 
fiber{1}(6).ROI = [21,37,34,5,18,29];

%DL102 180506
fiber{3}(1).ROI = [33,53,56,62,42,43,46,45,48,49]; 
fiber{3}(2).ROI = [32,28,36,30,38,37]; 
fiber{3}(3).ROI = [47,57,58]; 
fiber{3}(4).ROI = [8,12,16,21]; 
fiber{3}(5).ROI = [59,63,60,64]; 
fiber{3}(6).ROI = [10,19,6];

%DL102 180514
fiber{5}(1).ROI = [18,42,70,69,128,134,133,151,139,28,108,116,130]; 
fiber{5}(2).ROI = [42,78,105,90,91,99,115]; 
fiber{5}(3).ROI = [22,45,41,47,73,93,107,120]; 
fiber{5}(4).ROI = [103,20,26,39,21]; 
fiber{5}(5).ROI = [5,4,6,13,15]; 

%DL112_180528 - registration issues around CSD fixed october 2022

fiber{8}(1).ROI = [74,33,15,3,2]; %[148,63,31,46,35,76,38,48,95,112,118]; 
fiber{8}(2).ROI = [44,79]; %[148,137,87,128]; 
fiber{8}(3).ROI = [123,125];%[37,58,32]; 
fiber{8}(4).ROI = [65,43,36,34];%[41,42,39,52,68,93,104]; 
fiber{8}(5).ROI = [132,154,179,29,172,180,169,71]; %[167,212,186,207,174,178,149]; 
fiber{8}(6).ROI = [98,94,112,135,145];%[117,150,147,143,144];
fiber{8}(7).ROI = [25,47,9,30,31,48,66];%[98,110,105,133,175]; 
fiber{8}(8).ROI = [46,54,126,59];%[65,30,27,64,88,121,132,141]; 
fiber{8}(9).ROI = [26,28,75,80,100,115];%[109,103,115,97,82]; 
fiber{8}(10).ROI = [118,144,128,102,111,127,143];% [80,85]; 
fiber{8}(11).ROI = [84,134,136,165,164];% [123,146,131,160,129,153,155,158,154,122,169,168]; 
fiber{8}(12).ROI = [96,120,137,141,146,155];% [77,182,111,130,190]; 

% DL112 180617 FOV1
fiber{9}(1).ROI = [5,6,10,13,17,20,18,25,74,27,88,100];  % 
fiber{9}(2).ROI = [129,130]; 
fiber{9}(3).ROI = [164,137,128,107,83]; 
fiber{9}(4).ROI = [201,241,212,271,230]; 
fiber{9}(5).ROI = [103,123,165,152,198];
fiber{9}(6).ROI = [225,145,117,72,85,99,160,167,140,135,142,181]; 
fiber{9}(7).ROI = [257,261,202,182]; 
fiber{9}(8).ROI = [81,78,82,87,90,92,102,143,151]; 
fiber{9}(9).ROI = [97,182]; 
fiber{9}(10).ROI = [199,211,185,104,158,172,170,174,256]; 
fiber{9}(11).ROI = [122,274,260,254,227];  % 
fiber{9}(12).ROI = [285,284,302,290,319]; 
fiber{9}(13).ROI = [204,169]; 
fiber{9}(14).ROI = [11,21,35]; 
fiber{9}(15).ROI = [58,56];
fiber{9}(16).ROI = [16,136,110]; 
fiber{9}(17).ROI = [150,286]; 
fiber{9}(18).ROI = [31,28]; 
fiber{9}(19).ROI = [186,156,149]; 
fiber{9}(20).ROI = [216,138,208,245,134,141,154,220]; 
fiber{9}(21).ROI = [120,218,217];  % 
fiber{9}(22).ROI = [235,244,268,240]; 
fiber{9}(23).ROI = [265,252,250,213,109,80]; 
fiber{9}(24).ROI = [144,98]; 
fiber{9}(25).ROI = [90,92];
fiber{9}(26).ROI = [290,319]; 
fiber{9}(27).ROI = [1,19,42]; 
fiber{9}(28).ROI = [52,51,67,59]; 
fiber{9}(29).ROI = [47,121,184,263]; 
fiber{9}(30).ROI = [146,298,251,293,313]; 
fiber{9}(31).ROI = [36,96,116,307,50];  % 
fiber{9}(32).ROI = [54,93,205,68,26,37,299,247]; 
fiber{9}(33).ROI = [311,309,305,325,315]; 



%DL112_180617
%{
fiber{9}(1).ROI = [7,143,5,8,10,14,15,36,28,54,42];  % 143 54
fiber{9}(2).ROI = [5,85,79,47,65,55]; 
fiber{9}(3).ROI = [34,38,44,51,62,68,118,104]; 
fiber{9}(4).ROI = [197,164,137]; 
fiber{9}(5).ROI = [145,162,108,109,134,187];
fiber{9}(6).ROI = [110,102,86]; 
fiber{9}(7).ROI = [100,136,125,193,181]; 
fiber{9}(8).ROI = [112,127,198,217]; 
fiber{9}(9).ROI = [91,155,192,195,180]; 
fiber{9}(10).ROI = [13,17,21,27]; 
%}
%fiber{9}(5).ROI = [139,124,84,31,53,74,103,130,119,105,148]; % formerly fiber 5
%fiber{9}(8).ROI = [138,160,167,194,80,141,152,146,203]; % % formerly fiber 8

% DL115 180704
fiber{10}(1).ROI = [288,289,265,381,382]; 
fiber{10}(2).ROI = [147,148,151,157,160,161,162,163,168,172,174]; 
fiber{10}(3).ROI = [198,199,200,201]; 
fiber{10}(4).ROI = [208,211,213,215,217,218,219,221]; 
fiber{10}(5).ROI = [205,206,207,210,212,214,216]; 
fiber{10}(6).ROI = [237,240,242,245,247,248,249,252,254,258,260]; 
fiber{10}(7).ROI = [238,239,241,243,244,246,250,251,255,256]; 
fiber{10}(8).ROI = [113,114,121]; 
fiber{10}(9).ROI = [262,264]; 
fiber{10}(10).ROI = 356:362; 
fiber{10}(11).ROI = [367,368,371,373,374,375,376]; 
fiber{10}(12).ROI = [377:380]; 
fiber{10}(13).ROI = [446:450]; 
fiber{10}(14).ROI = [400:403,406,407]; 
fiber{10}(15).ROI = [410:419]; 
fiber{10}(16).ROI = [420,421]; 
fiber{10}(17).ROI = [422:426]; 
fiber{10}(18).ROI = [477,479,480]; 
fiber{10}(19).ROI = [482:486]; 


% DL115 180707 - missing frames around CSD make this data problematic for CSD fixed 10/22
% {
fiber{11}(1).ROI = [187,164,181,170,178,173,176]; %[43,15,16,25,46,32,60]; % [43,60];
fiber{11}(2).ROI = [219,242,256]; %[566,626,599,405,395,504,483,468,586,644,491]; % ,463 ,433 ,575 ,400
fiber{11}(3).ROI = [3,13,15,20]; %[7,30,39,9,12,6]; 
fiber{11}(4).ROI = [66,158]; %[47,42,84,156,176]; 
fiber{11}(5).ROI = [179,251,267,272,259,258]; %[417,501,472,662,636,625,412,326]; 
fiber{11}(6).ROI = [91,47,21,6,35]; %[460,525,385,389,349]; 
fiber{11}(7).ROI = [95,34]; %[512,377,351,141,140]; 
fiber{11}(8).ROI = [18,16];
fiber{11}(9).ROI = [24,8,11,27];
fiber{11}(10).ROI = [236,235];
fiber{11}(11).ROI = [38,55,29,25];
fiber{11}(12).ROI = [56,46,54];
fiber{11}(13).ROI = [238,194];
fiber{11}(14).ROI = [260,217,148];
fiber{11}(15).ROI = [145,264,209];
fiber{11}(16).ROI = [169,165,161];
%}
% DL117 170704 - no CSD
fiber{12}(1).ROI = [9, 2, 1, 12];
fiber{12}(2).ROI = [18, 14, 6, 7];
fiber{12}(3).ROI = [3, 4];

% DL118 180628 - missing key onset scans
%{
fiber{13}(1).ROI = [271,286,290,273,285,272,278,254,239,256]; % ,253 251,194
fiber{13}(2).ROI = [136,96,54,52,43,23];
fiber{13}(3).ROI = [73,170,226,208,188,162,142];
fiber{13}(4).ROI = [152,140,88,107,89,80,64,49];
fiber{13}(5).ROI = [292,287,211,201,198,160,144,105,130,131];
fiber{13}(6).ROI = [300,302,297];
fiber{13}(7).ROI = [134,94];
%}

% DL118 180715 - CSD is bad - fixed 10/22
fiber{14}(1).ROI = [23,35,53,82]; %[9,48,56,6,41,47,43,31,4]; 
fiber{14}(2).ROI = [76,48,55,80]; %[106,74,87,92,96];
fiber{14}(3).ROI = [118,112]; %[30,18,7];
fiber{14}(4).ROI = [102,99]; %[105,97,52,12];
fiber{14}(5).ROI = [122,107,115]; %[130,93,67];
fiber{14}(6).ROI = [63,84,78,19,27,20,14]; %[39,63,78,98,115,140];
fiber{14}(7).ROI = [60,33,16]; %[72,63];
fiber{14}(8).ROI = [47,21,31]; %[130,114,85];
fiber{14}(9).ROI = [93,88,9,117,105,116,109,110,95,108]; %[91,49,32];
fiber{14}(10).ROI = [111,15,41,46]; % [111,127,131,135,109];
fiber{14}(11).ROI = [103,57,113,28]; %
fiber{14}(12).ROI = [36,52,13,103]; %
fiber{14}(13).ROI = [74,83,92,75]; %
fiber{14}(14).ROI = [74,83,92,75]; %
fiber{14}(15).ROI = [50,49]; %

% DL67 170323 FOV2
%fiber{22}(1).ROI = [];

% DL122 180805
fiber{16}(1).ROI = [646,1177,1168,526];
fiber{16}(2).ROI = [506,1099,658,1250,635,644,1257,1174,1107,585,925,1106,670];
fiber{16}(3).ROI = [544,305,256];
fiber{16}(4).ROI = [416,894,452,297];
fiber{16}(5).ROI = [197,211,253,474,746,967,684,594];
fiber{16}(6).ROI = [22,113,93,71,202,54,29];
fiber{16}(7).ROI = [210,665];
fiber{16}(8).ROI = [237,259,314];
fiber{16}(9).ROI = [1028,835,1189,1349,1371,819,914,633,601,933,521,603,898];
fiber{16}(10).ROI = [691,1001,983];
fiber{16}(11).ROI = [1227,1118,1220,1201,1073];
fiber{16}(12).ROI = [1357,1322];
fiber{16}(13).ROI = [945,728,736,681,530,522,371,1086,735];
fiber{16}(14).ROI = [457,500];
fiber{16}(15).ROI = [100,44];
fiber{16}(16).ROI = [332,339,565];
fiber{16}(17).ROI = [1199,1077,1094];
fiber{16}(18).ROI = [380,352];
fiber{16}(19).ROI = [325];
fiber{16}(20).ROI = [1284,1302,1275,1320,1341,1251,1308];
fiber{16}(21).ROI = [536,598,969,1146,1314,1301];
fiber{16}(22).ROI = [1340,242];
fiber{16}(23).ROI = [1193,772,1133,751];
fiber{16}(24).ROI = [324,462,258];
fiber{16}(25).ROI = [393,373, 739,707,1115,1216];
fiber{16}(26).ROI = [367,593,563,566,568,628,699,1195,1187,529];
fiber{16}(27).ROI = [213,209,294];
fiber{16}(28).ROI = [1153,1152];
fiber{16}(29).ROI = [1192,1219,1317,1295];
fiber{16}(30).ROI = [1305,1288,1252,719,531,264,1151,287,226];
fiber{16}(31).ROI = [341,370,638,730,1212,842];
fiber{16}(32).ROI = [1288,1246];
fiber{16}(33).ROI = [47,40,36,35];
fiber{16}(34).ROI = [1345,1363,1358,1347];
fiber{16}(35).ROI = [676,1149,1148];
fiber{16}(36).ROI = [909,972,940,868,857,870,861];
fiber{16}(37).ROI = [331,636,710];
fiber{16}(38).ROI = [1138,1161,838,997];
fiber{16}(39).ROI = [282];
fiber{16}(40).ROI = [1229,1268,1315,1316,1318];
fiber{16}(41).ROI = [881,401,629,491,363,518,904,740,1034,654,995,558,931,738,1244,1213,655,1114,1155];
fiber{16}(42).ROI = [503,346,366,442,446,903,575,978,450,1116];
fiber{16}(43).ROI = [348,291,260,243,317,321];
fiber{16}(44).ROI = [1366,711,1160];
fiber{16}(45).ROI = [479,555,501,609];
fiber{16}(46).ROI = [960,1046];
fiber{16}(47).ROI = [790,157,177,185,246,41,14,25];
fiber{16}(48).ROI = [351,229,238];
fiber{16}(49).ROI = [1287,519,543];
fiber{16}(50).ROI = [343,407,319];
fiber{16}(51).ROI = [556];
fiber{16}(52).ROI = [579];
fiber{16}(53).ROI = [292,1045];
fiber{16}(54).ROI = [1194,1173,1249,1214,1267];
fiber{16}(55).ROI = [741,276,240];
fiber{16}(56).ROI = [1175,1178];
fiber{16}(57).ROI = [752,766,388];
fiber{16}(58).ROI = [1277,1209,586,1164,1081,649];
fiber{16}(59).ROI = [159,122];
fiber{16}(60).ROI = [218,263,398];
fiber{16}(61).ROI = [1292,1247,1090,1280,987];
fiber{16}(62).ROI = [577,512,441,996,408,303,106,231,245,1018];
fiber{16}(63).ROI = [743,1232,1270,1218,727,1191,1190,428,1279,1196,1343,1296,1321,657,1256,1109,770,417];
fiber{16}(64).ROI = [42,32,26];
fiber{16}(65).ROI = [716,731,806];
fiber{16}(66).ROI = [255,108,73,327,272,270];
fiber{16}(67).ROI = [158,131,98,85,60,89,94,109,195,195,117,233,230,167];
fiber{16}(68).ROI = [461,1035,1265,1323];
fiber{16}(69).ROI = [313,1065,391,959];
fiber{16}(70).ROI = [1089,708,1206,1228,470,1171];



%
fiber{27}(1).ROI = [73,120,77];
fiber{27}(2).ROI = [49,50];
fiber{27}(3).ROI = [25,26];
fiber{27}(4).ROI = [4,5,6];
fiber{27}(5).ROI = [17:19];
fiber{27}(6).ROI = [17,63,64];
fiber{27}(7).ROI = [129,38];
fiber{27}(8).ROI = [68:70];
fiber{27}(9).ROI = [60,61,106,29];
fiber{27}(10).ROI = [55,56];
fiber{27}(11).ROI = 42:43;
fiber{27}(12).ROI = [54,122,119];
fiber{27}(13).ROI = [11:14];
fiber{27}(14).ROI = [132,117];
fiber{27}(15).ROI = [110,138];
fiber{27}(16).ROI = [30:32];
fiber{27}(17).ROI = [78:79];
fiber{27}(18).ROI = [7:10];
fiber{27}(19).ROI = [123,51,52];
fiber{27}(20).ROI = [95,101];
fiber{27}(21).ROI = [36,37];
fiber{27}(22).ROI = [22:24,82,71,72]; % 168,
fiber{27}(23).ROI = [27,20]; % 21
fiber{27}(24).ROI = [47,48];
fiber{27}(25).ROI = [15,16];

% DL67 170601
fiber{28}(1).ROI = [79:83];
fiber{28}(2).ROI = [119:122];
fiber{28}(3).ROI = [62:64];
fiber{28}(4).ROI = [98:101];
fiber{28}(5).ROI = [369:371];
fiber{28}(6).ROI = [222:225];
fiber{28}(7).ROI = [215:217];
fiber{28}(8).ROI = [395:397,109,110,5,65];
fiber{28}(9).ROI = [434:435];
fiber{28}(10).ROI = [346:347];
fiber{28}(11).ROI = 87:88;

% DL68 170523
fiber{30}(1).ROI = [46,54,67,81,78,84,82];
fiber{30}(2).ROI = [29,24,18,35];
fiber{30}(3).ROI = [15,31,53,72,80,56,70];
fiber{30}(4).ROI = [108,103,102,95,83,58];
fiber{30}(5).ROI = [49,38,17,6];
fiber{30}(6).ROI = [30,34,12,25];
fiber{30}(7).ROI = [40,37,61];
fiber{30}(8).ROI = [77,43,66,71];
fiber{30}(9).ROI = [91,98];
fiber{30}(10).ROI = [32,48,47];

% DL68 170712
fiber{31}(1).ROI = [29,35,37,8,30,34,26,27]; %[35,44,56,75,71];
fiber{31}(2).ROI = [5,23,10]; %[67,52,71,70,63,58];
fiber{31}(3).ROI = [17,7,3]; %[58,64,65,59];
fiber{31}(4).ROI = [33,19]; %[35,72,74,50];
fiber{31}(5).ROI = [12,31,28,14,9,16,21,1]; %[25,20,9,12];
fiber{31}(6).ROI = [24,36,20]; %[18,13,10,15];


% DL72 170429
fiber{33}(1).ROI = [278,279,372,273,373,280,491,492,281,282,420,283]; % [278, 279,280,491,492,281,282,420,283];  % ,374
fiber{33}(2).ROI = [131,132,133];
fiber{33}(3).ROI = [59,60,61,117];
fiber{33}(4).ROI = [254,255,256,257,258];
fiber{33}(5).ROI = [512,513];
fiber{33}(6).ROI = [21,23];
fiber{33}(7).ROI = [29,30,31,33];
fiber{33}(8).ROI = [469,79,478,479,48];
fiber{33}(9).ROI = [413,412,414,416];
fiber{33}(10).ROI = [36,37];
fiber{33}(11).ROI = [120,121,122];
fiber{33}(12).ROI = [266,267,269,270];
fiber{33}(13).ROI = [110,111,113,114,116];
fiber{33}(14).ROI = [284,286,287,289];
fiber{33}(15).ROI = [212:215];
fiber{33}(16).ROI = [215:217];
fiber{33}(17).ROI = [85:90];
fiber{33}(18).ROI = [24:26];
fiber{33}(19).ROI = [160,161,164];
fiber{33}(20).ROI = [17:20,22];
fiber{33}(21).ROI = [175,177];
fiber{33}(22).ROI = [174,176,179];
fiber{33}(23).ROI = [343:345,347];
fiber{33}(24).ROI = [153,156,157,158];
fiber{33}(25).ROI = [178,181,182,180];
fiber{33}(26).ROI = [398:401];
fiber{33}(27).ROI = [432,433];
fiber{33}(28).ROI = [383:386];
fiber{33}(29).ROI = [372:375];
fiber{33}(30).ROI = [472,473,105];
fiber{33}(31).ROI = [119,123,124];

% DL72 170609
fiber{35}(1).ROI = [389:392];
fiber{35}(2).ROI = [155,157,160,162,165];
fiber{35}(3).ROI = [172:176];
fiber{35}(4).ROI = [45:48,52,57];
fiber{35}(5).ROI = [58,53:56,59];
fiber{35}(6).ROI = [102,106,110,112];
fiber{35}(7).ROI = [100,102,103,104,105];
fiber{35}(8).ROI = [201,203:208,213];
fiber{35}(9).ROI = [238:242];
fiber{35}(10).ROI = [243:251];
fiber{35}(11).ROI = [74:77];
fiber{35}(12).ROI = [131:135,140];
fiber{35}(13).ROI = [138:139,143,144];
fiber{35}(14).ROI = [287:290];
fiber{35}(15).ROI = [60:63,65,70,71,72];
fiber{35}(16).ROI = [417:419,421:425,427];
fiber{35}(17).ROI = [87:91,94];
fiber{35}(18).ROI = [92,93,95:98];
fiber{35}(19).ROI = [266,270:274];
fiber{35}(20).ROI = [306:310];
fiber{35}(21).ROI = [321,323,324,326,327];
fiber{35}(22).ROI = [320,322,325,328,330,331];
fiber{35}(23).ROI = [346:349,351,352];
fiber{35}(24).ROI = [216:224];
fiber{35}(25).ROI = [254,255,257,260,261,262];
fiber{35}(26).ROI = [275,276,278:284];
fiber{35}(27).ROI = [336,338:340];
fiber{35}(28).ROI = [311,312,315,316,318];
fiber{35}(29).ROI = [371,373,375];
fiber{35}(30).ROI = [408:412];
fiber{35}(31).ROI = [413,414];
fiber{35}(32).ROI = [457:459];
fiber{35}(33).ROI = [400:406];
fiber{35}(34).ROI = [491:497];
fiber{35}(35).ROI = [505:510];
fiber{35}(36).ROI = [433:436];
fiber{35}(37).ROI = [465:469];
fiber{35}(38).ROI = [470:474];
fiber{35}(39).ROI = [526:530];
fiber{35}(40).ROI = [389:392];
fiber{35}(41).ROI = [428:430];
fiber{35}(42).ROI = [511:514];
fiber{35}(43).ROI = [502:504];
fiber{35}(44).ROI = [533:535];

% DL72 170614
fiber{36}(1).ROI = [6:10];
fiber{36}(2).ROI = [304,324:326];
fiber{36}(3).ROI = [397:402];
fiber{36}(4).ROI = [442:443,455];
fiber{36}(5).ROI = [492:493];
fiber{36}(6).ROI = [245,246,252];
fiber{36}(7).ROI = [368:370];
fiber{36}(8).ROI = [305,377:379];
fiber{36}(9).ROI = [342:345];
fiber{36}(10).ROI = [244:246,251:252];
fiber{36}(11).ROI = [12:13];
fiber{36}(12).ROI = [15:16];
fiber{36}(13).ROI = [41,42];
fiber{36}(14).ROI = [409,410];
fiber{36}(15).ROI = [70,72:74];
fiber{36}(16).ROI = [255:256];
fiber{36}(17).ROI = [459:463];
fiber{36}(18).ROI = [413:417];
fiber{36}(19).ROI = [444:451];
fiber{36}(20).ROI = [63,65:68];
fiber{36}(21).ROI = [157:166];
fiber{36}(22).ROI = [304:306];
fiber{36}(23).ROI = [167:171];
fiber{36}(24).ROI = [287,289:293,296];
fiber{36}(25).ROI = [92,94:95,93,98,97,96];
fiber{36}(26).ROI = [92,88,89,90,91];
fiber{36}(27).ROI = [113,114,116:119];
fiber{36}(28).ROI = [131:134];
fiber{36}(29).ROI = [272:275];
fiber{36}(30).ROI = [196:198,200];
fiber{36}(31).ROI = [79:82,84];
fiber{36}(32).ROI = [183:188];
fiber{36}(33).ROI = [315:323];
fiber{36}(34).ROI = [380:384];
fiber{36}(35).ROI = [517:518];

% DL75
fiber{37}(1).ROI = [17,23,42,44,36,54,59,67,65,66,70,76,72,75,73,80,74,84]; 
fiber{37}(2).ROI = [18,15,20,25,19,14,12,35]; 
fiber{37}(3).ROI = [83,57,40,1]; 
fiber{37}(4).ROI = [77,51,81,69,61]; 
fiber{37}(5).ROI = [33,30,34,37,45,39,49,60]; 
fiber{37}(6).ROI = [7,4,3,2,5,58]; 
fiber{37}(7).ROI = [68,62,33,56,55,41,26]; 
fiber{37}(8).ROI = [24,28,27]; 

% DL89 171113
fiber{38}(1).ROI = [333,281,260,263];
fiber{38}(2).ROI = [92,104,121,99];
fiber{38}(3).ROI = [2,5,28,25,62];
fiber{38}(4).ROI = [131,106];
fiber{38}(5).ROI = [76,134,255];
fiber{38}(6).ROI = [256,224,236,145,82,187,296];
fiber{38}(7).ROI = [348,403,408,392,385,397,387,294];
fiber{38}(8).ROI = [401,411,379,313];
fiber{38}(9).ROI = [87,124,219,161,169,226,225,264];
fiber{38}(10).ROI = [94,81,70,169];
fiber{38}(11).ROI = [63,100,203,189];
fiber{38}(12).ROI = [152,151,97,49,24];
fiber{38}(13).ROI = [50,69,37,44,40,49];
fiber{38}(14).ROI = [278,346,352,231,275,271];
fiber{38}(15).ROI = [295,306,292,300,245,142,159,174,208];
fiber{38}(16).ROI = [407,418,386,363,372,380,395,421,406];
fiber{38}(17).ROI = [23,21,19,26,27,36,30,20];
fiber{38}(18).ROI = [35,57,65,74,75];
fiber{38}(19).ROI = [149,153,154,146,209];
fiber{38}(20).ROI = [378,336,362,358];

% DL89 171116 
fiber{39}(1).ROI = [62, 1, 2, 5, 12]; 
fiber{39}(2).ROI = [7, 8];
fiber{39}(3).ROI = [10, 11]; 
fiber{39}(4).ROI = [31, 28];
fiber{39}(5).ROI = [24, 25];
fiber{39}(6).ROI = [14, 15];
fiber{39}(7).ROI = [49, 79];

% DL89 171122
fiber{41}(1).ROI = [96,99,65,62,57,37,38,41,52]; 
fiber{41}(2).ROI = [36,30,14];
fiber{41}(3).ROI = [102,112,134,80,66,50];
fiber{41}(4).ROI = [192,183,182,170,82,117,121]; %
fiber{41}(5).ROI = [140,114,103]; 
fiber{41}(6).ROI = [91,150,58,87,118,95,72,79,59];
fiber{41}(7).ROI = [49,35,31,18,10,3];
fiber{41}(8).ROI = [129,169,151,83,161];
fiber{41}(9).ROI = [90,86,63,56];
fiber{41}(10).ROI = [138,162,178,177,168,109];
fiber{41}(11).ROI = [144,141,155,158,148,147];
fiber{41}(12).ROI = [145,156,173,107];
fiber{41}(13).ROI = [22,33,28];

% DL89_20171119
fiber{43}(1).ROI = [34:37]; 
fiber{43}(2).ROI = [106,109,558,559,110,111];
fiber{43}(3).ROI = [89,92,93,96];
fiber{43}(4).ROI = [84,85,90]; %
fiber{43}(5).ROI = [139:141]; 
fiber{43}(6).ROI = [132,133,135:138];
fiber{43}(7).ROI = [191,193];
fiber{43}(8).ROI = [706:708];
fiber{43}(9).ROI = [719:721];
fiber{43}(10).ROI = [470:472];
fiber{43}(11).ROI = [476:478];
fiber{43}(12).ROI = [267,266,268,269];
fiber{43}(13).ROI = [570:572];
fiber{43}(14).ROI = [210:213];
fiber{43}(15).ROI = [217:220];
fiber{43}(16).ROI = [222:225];
fiber{43}(17).ROI = [779:780];
fiber{43}(18).ROI = [50,51,56];
fiber{43}(19).ROI = [179:182];
fiber{43}(20).ROI = [53,52,54,56,57];
fiber{43}(21).ROI = [678,241,418:420];

Nfiber = cellfun(@numel, fiber);
xFiber = find(Nfiber>0);
for x = intersect(xPresent, xFiber)
    fiber{x} = MakeFibers(expt{x}, fiber{x}, ROI{x}, false); %true false
    expt{x}.sbx.interp = strcat(expt{x}.dir, expt{x}.name, '.sbx_interp ');
    fluor{x} = GetAxonFluor(expt{x}, catInfo{x}, fiber{x}, fluor{x}, deform{x}, 'window',find(Tscan{x}{1}<=32,1,'last'), 'overwrite',false);
end

%% Correlation vs distance: intra- vs extra-fiber ROIs 
intraTable = cell(1,Nexpt); extraTable = cell(1,Nexpt); % [];

sepBinWidth = 20;
sepEdges = 0:sepBinWidth:700; Nsep = numel(sepEdges);
sepCorrBins = cell(Nsep,2);
for x = intersect(xFiber, xPresent)
    if ~isnan(expt(x).csd)
        fluorCat = [fluor{x}(1:expt(x).csd-1).dFF];
        fluorCat = vertcat(fluorCat.ROI);
    else
        fluorCat = [fluor{x}.dFF];
        fluorCat = vertcat(fluorCat.ROI);
    end
    % Pairwise distance
    tempCent = vertcat(ROI{x}.cent);
    roiSep = expt(x).umPerPixel*squareform(pdist(tempCent(:,1:2)));
    %imagesc(roiSep); colorbar; axis square;
    
    % Pairwise correlation
    fluorCatCorr = corr(fluorCat, 'Rows','complete');
    %imagesc(fluorCatCorr); colorbar; axis square;
    
    for f = 1:Nfiber(x)
        tempExtraROI = setdiff(1:expt(x).Nroi, fiber{x}(f).ROI);
        for roi = fiber{x}(f).ROI
            tempIntraROI = setdiff(fiber{x}(f).ROI, roi);
             
            intraTable{x} = vertcat( intraTable{x}, [roiSep(tempIntraROI,roi), fluorCatCorr(tempIntraROI,roi)] );
            extraTable{x} = vertcat( extraTable{x}, [roiSep(tempExtraROI,roi), fluorCatCorr(tempExtraROI,roi)] );
        end
        %{
        tempIntraPair = nchoosek(fiber{x}(f).ROI, 2);
        for p = 1:size(tempIntraPair,1)
            k = k+1;
            intraTable(k,1) = roiSep(tempIntraPair(p,1), tempIntraPair(p,2));
            intraTable(k,2) = fluorCatCorr(tempIntraPair(p,1), tempIntraPair(p,2));
            %plot(roiSep(tempIntraPair(p,1), tempIntraPair(p,2)), fluorCatCorr(tempIntraPair(p,1), tempIntraPair(p,2)), '.'); hold on;
        end
        
        for roi = fiber{x}(f).ROI
            extraTable(:,1) = roiSep(roi, setdiff(1:expt(x).Nroi, fiber{x}(f).ROI));
        end
        %}
    end
    %{
    plot(extraTable{x}(:,1), extraTable{x}(:,2), 'r.' ); hold on;
    plot(intraTable{x}(:,1), intraTable{x}(:,2), 'k.', 'MarkerSize',10 ); 
    ylabel('Pairwise Correlation'); xlabel('Pairwise Distance (um)');
    pause;
    %}
    %tempTable = [intraTable{x}(1:10,:), discretize(intraTable{x}(1:10,1), sepEdges)];
    %[sepBinInd, ~, sepBinPair] = unique(discretize(intraTable{x}(1:10,1), sepEdges));
    %histcounts(intraTable{x}(:,1), sepEdges);
    %histogram(intraTable{x}(:,1), sepEdges, 'Normalization','probability');
    
    %{
    intraSepBins = discretize(intraTable{x}(:,1), sepEdges);
    extraSepBins = discretize(extraTable{x}(:,1), sepEdges);
    for bin = 1:max(intraSepBins)
        tempIntraBinROI = find(intraSepBins == bin)';
        tempExtraBinROI = find(extraSepBins == bin)';
        sepCorrBins{bin,1} = intraTable{x}(tempIntraBinROI,2);
        sepCorrBins{bin,2} = extraTable{x}(tempExtraBinROI,2);
    end
    %plot( cellfun(@mean, sepCorrBins) )
    %}
    
end

intraTablePool = vertcat( intraTable{xPresent} );
extraTablePool = vertcat( extraTable{xPresent} );
intraSepBins = discretize(intraTablePool(:,1), sepEdges);
extraSepBins = discretize(extraTablePool(:,1), sepEdges);
for bin = 1:max(intraSepBins)
    tempIntraBinROI = find(intraSepBins == bin)';
    tempExtraBinROI = find(extraSepBins == bin)';
    sepCorrBins{bin,1} = intraTablePool(tempIntraBinROI,2);
    sepCorrBins{bin,2} = extraTablePool(tempExtraBinROI,2);
end
sepCorrBinN = cellfun(@numel, sepCorrBins);
sepCorrBins = sepCorrBins(sepCorrBinN(:,1) > 0,:); % NaNs mess up plotshaded
sepCorrBinN = cellfun(@numel, sepCorrBins);
sepCorrBinMean = cellfun(@mean, sepCorrBins);
sepCorrBinStd = cellfun(@std, sepCorrBins);
sepBinCent = (sepBinWidth/2)*(1:size(sepCorrBins,1));
sepCorrBinSEM = sepCorrBinStd./sqrt(sepCorrBinN);
%plot(1:size(sepCorrBins,1), sepCorrBinMean);


FiberCorrVsSeparation = figure('Units','inches', 'Position', [6, 5, 1.1, 1.1], 'Color','w'); 
sepCorrBinError = sepCorrBinSEM; sepCorrBinStd; %
plotshaded(sepBinCent, [sepCorrBinMean(:,1)-sepCorrBinError(:,1), sepCorrBinMean(:,1), sepCorrBinMean(:,1)+sepCorrBinError(:,1)]', '-'  ); hold on;
plotshaded(sepBinCent, [sepCorrBinMean(:,2)-sepCorrBinError(:,2), sepCorrBinMean(:,2), sepCorrBinMean(:,2)+sepCorrBinError(:,2)]', '-', 'r'  );
set(gca, 'Xtick',0:50:250, 'Yscale','log', 'TickDir','out', 'box','off', 'Position',[0.02,0.02,0.96,0.96], 'TickLength',[TL,0]); %,  , 'Ytick',[0.1,1]
%ylim([0,0.6]);
%xlabel('Separation (\mum)', 'Interpreter','tex');
%ylabel('Correlation');
figPath = sprintf('%sFiberCorrVsSeparation.tif', fig4Dir);
exportgraphics(FiberCorrVsSeparation, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);

%%
for x = 10%37
    for run = 1 %3
        for f = 1
            plot(Tscan{x}{run}, fluor{x}(run).z.ROI(:,fiber{x}(f).ROI)); hold on;
            for e = find([stillEpoch{x}.run] == run)
                line( [stillEpoch{x}(e).Tstart, stillEpoch{x}(e).Tstop], -2*[1,1], 'color','r'); 
            end
            axis tight;
            pause;
        end
    end
end

%%
% Search for discordant events within a single fiber
FS = 8;
close all
figure('WindowState','maximized');
for x = 33%intersect(xPresent, xFiber)
    for f = 1 %1:Nfiber(x)
        SP(1) = subplot(2,2,1); cla;
        imshow( expt(x).maxProj, prctile(expt(x).maxProj(:), [5,99.5]) ); hold on;
        for roi = fiber{x}(f).ROI
            plot( ROI{x}(roi).footprintEdge(:,2), ROI{x}(roi).footprintEdge(:,1), 'b.', 'MarkerSize',2 ) 
        end
        
        SP(2) = subplot(2,2,3); cla;
        imshow( label2rgb(fiber{x}(f).labelFoot ), [] ); %hold on;
        for roi = fiber{x}(f).ROI
            text( ROI{x}(roi).cent(1), ROI{x}(roi).cent(2), num2str(roi), 'color','m', 'FontWeight','bold', 'HorizontalAlignment','center', 'FontSize',FS )
        end
        impixelinfo;
        linkaxes(SP,'xy');
        
        for run = 1 %2
            Froi = {fluor{x}(run).z.ROI(:,flip(fiber{x}(f).ROI))};
            Froi_spread = CalSpread( Tscan{x}(run), Froi, 'show',false, 'SepPct',85 );
            Fmean = mean(Froi{1}, 2, 'omitnan');
            Fsub = Froi{1} - Fmean;
            Fsub_spread = CalSpread( Tscan{x}(run), {Fsub}, 'show',false, 'SepPct',85 ); % fluor{x}(run).z.ROI(:,fiber{x}(f).ROI)
            
            sp(1) = subplot(4,2,2); cla;
            plot( Froi_spread{1} );
            ylim([-Inf, Inf]);
            title( sprintf('[x,f,run] = [%i, %i, %i]' , x, f, run) );
            ylabel('ROI traces');
            
            sp(2) = subplot(4,2,4); cla;
            plot( Fsub_spread{1} );
            ylim([-Inf, Inf]);
            ylabel('Mean-subtracted traces');
            
            sp(3) = subplot(4,2,6); cla;
            imagesc(eventRaster{x}{run}(:,fiber{x}(f).ROI)')
            set(gca, 'Ytick',1:fiber{x}(f).Nroi, 'YtickLabel',fiber{x}(f).ROI, 'XTickLabel',[]) % 
            title('Events');
            
            sp(4) = subplot(4,2,8); cla;
            plot( sum(eventRaster{x}{run}(:,fiber{x}(f).ROI), 2) );
            ylabel('# of ROI');

            linkaxes(sp, 'x');
            axis tight;
            
            pause;
        end
    end
end

%%  Write movies of fibers during periods of discordance

fiberDir = 'C:\Users\ablaeser\Documents\Afferent Paper\PaperFigures\Movies\';

X = 33;
F = 1; 
R = 1;
S = expt(X).scanLims(R) + 30000:31500;


%{
X = 28;
F = 2; 
R = 1;
S = expt(X).scanLims(R) + 25001:26500;
%}
cropPad = [10,10];
discordName = sprintf('%s_fiber%i_run%i_scans%i-%i', expt(X).name, F, R, S(1), S(end));
fiberBoxX = vertcat(ROI{X}(fiber{X}(F).ROI).box_x);
fiberBoxY = vertcat(ROI{X}(fiber{X}(F).ROI).box_y);
fiberCrop = round([min(fiberBoxX(:,1)), expt(X).Ncol - max(fiberBoxX(:,2)), min(fiberBoxY(:,1)), expt(X).Nrow - max(fiberBoxY(:,2))] + [-cropPad(1), cropPad(1), -cropPad(2), cropPad(2)]);
%round([min(fiberBoxX(:,1)), max(fiberBoxX(:,2)), min(fiberBoxY(:,1)), max(fiberBoxY(:,2))] + [-cropPad(1), cropPad(1), -cropPad(2), cropPad(2)]);

WriteSbxPlaneTif(expt(X).sbx, catInfo(X), 1, 'dir',fiberDir, 'name',discordName, 'type','aff', 'edge',fiberCrop, 'firstScan',S(1), 'Nscan',numel(S), 'overwrite',true ); % segParams{X}.edges


%% Show discordant period within a single fiber
FS = 8;
close all
figure('WindowState','maximized');
SP(1) = subplot(2,2,1); cla;
imshow( expt(X).maxProj, prctile(expt(X).maxProj(:), [5,99.5]) ); hold on;
for roi = fiber{X}(F).ROI
    plot( ROI{X}(roi).footprintEdge(:,2), ROI{X}(roi).footprintEdge(:,1), 'b.', 'MarkerSize',2 )
end

SP(2) = subplot(2,2,3); cla;
imshow( label2rgb(fiber{X}(F).labelFoot ), [] ); %hold on;
for roi = fiber{X}(F).ROI
    text( ROI{X}(roi).cent(1), ROI{X}(roi).cent(2), num2str(roi), 'color','m', 'FontWeight','bold', 'HorizontalAlignment','center', 'FontSize',FS )
end
impixelinfo;
linkaxes(SP,'xy');
%xlim(fiberCrop(1:2)); ylim(fiberCrop(3:4));

Froi = {fluor{X}(R).z.ROI(S,flip(fiber{X}(F).ROI))};
Froi_spread = CalSpread({1:numel(S)}, Froi, 'show',false, 'SepPct',85 ); % Tscan{X}(R)
Fmean = mean(Froi{1}, 2, 'omitnan');
Fsub = Froi{1} - Fmean;
Fsub_spread = CalSpread( {1:numel(S)}, {Fsub}, 'show',false, 'SepPct',85 ); % fluor{X}(R).z.ROI(:,fiber{X}(F).ROI)

sp(1) = subplot(4,2,2); cla;
plot( Froi_spread{1} );
ylim([-Inf, Inf]);
title( sprintf('[X,F,R] = [%i, %i, %i]' , X, F, R) );
ylabel('ROI traces');

sp(2) = subplot(4,2,4); cla;
plot( Fsub_spread{1} );
ylim([-Inf, Inf]);
ylabel('Mean-subtracted traces');

sp(3) = subplot(4,2,6); cla;
imagesc(eventRaster{X}{R}(S,fiber{X}(F).ROI)')
set(gca, 'Ytick',1:fiber{X}(F).Nroi, 'YtickLabel',fiber{X}(F).ROI, 'XTickLabel',[]) %
title('Events');

sp(4) = subplot(4,2,8); cla;
plot( sum(eventRaster{X}{R}(S,fiber{X}(F).ROI), 2) );
ylabel('# of ROI');

linkaxes(sp, 'x');
axis tight;

%% Plot fiber level fluor vs individual ROIs
close all;
for x = 5 % xFiber
    Tcat = vertcat(Tscan{x}{:});
    tempFluor = [fluor{x}.dFF];
    catFluor_ROI = vertcat(tempFluor.ROI);
    catFluor_fiber = vertcat(tempFluor.axon);
    figure;
    for f = 1:Nfiber(x)
        clf;
        subplot(1,2,1);
        imshow(label2rgb(fiber{x}(f).labelFoot)) 
        subplot(1,2,2)
        plot(Tcat, catFluor_ROI(:,fiber{x}(f).ROI), 'color',[1,0,0,0.2] ); hold on;
        plot(Tcat, catFluor_fiber(:,f), 'color',[0,0,0,0.7], 'lineWidth',1.5 );
        %fiber{x}(f).ROI
        pause;
    end
end