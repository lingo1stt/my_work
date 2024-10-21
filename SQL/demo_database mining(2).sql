#test
show databases;
#使用該資料級
use field_trial;
show tables;
##描述資料結構
describe Seed_Info ; 

SELECT * FROM Company_Info;
SELECT Parent_company, Parent_company  from Company_Info ci ;
SELECT * FROM  Seed_Info si  WHERE  Maturity  < 131;
SELECT * FROM  Seed_Info si  
WHERE  Hybrid_name = "P0720Q";
###選取多個name
SELECT * FROM  Seed_Info si  
WHERE  Hybrid_name IN ("P0720Q","DKC59-82RIB");

## %代表任何長度的東西
SELECT * FROM Seed_Info si  WHERE  Hybrid_name  LIKE "P%1%";
## _代表長度為1的任何東西
SELECT * FROM Seed_Info si  WHERE  Hybrid_name  LIKE "P_1%";

#############################################
##將SEED INFO 以及Company_info兩表格的column合併，而row之間的關係由company iD來確認
SELECT * FROM Seed_Info  JOIN Company_Info  
ON Seed_Info.Company_id = Company_Info.Company_ID 
WHERE Seed_Info.Maturity > 110 and Seed_Info.Maturity <113;
##記數(nrow)
SELECT count(*) from Company_Info join Seed_Info ;

###選取不重複的row
SELECT DISTINCT Company_Info.
Parent_company
FROM Seed_Info JOIN Company_Info
ON Seed_Info.Company_id =
Company_Info.Company_ID;

##limit:只取前n行
SELECT seed_id, maturity from Seed_Info limit 5;

###order:排序
SELECT * FROM  Seed_Info  order by Maturity LIMIT 5;
SELECT * FROM  Seed_Info  order by Company_id  DESC  limit 10;
SELECT * FROM  Seed_Info  order by Company_id  ASC  limit 10;

###多重排序(先根據companyid排序，在同一company id中又根據Maturity排序)
SELECT * FROM  Seed_Info  order by Company_id  DESC,maturity asc

###GROUP BY

SELECT State_Info.State_name
FROM State_Info JOIN Farm_Info
ON State_Info.State_ID = Farm_Info.State_ID
GROUP BY State_Info.State_name, Farm_Info.Irrigation  ;
###having
SELECT State_Info.State_name
FROM State_Info JOIN Farm_Info
ON State_Info.State_ID = Farm_Info.State_ID
GROUP BY State_Info.State_name, Farm_Info.Irrigation 
HAVING count(*) >=2;
###根據state name出現的次數進行統計(如果加上FARM INFO.irrigation則代表兩指標都要相同才會被分在同一群)
#select * 代表可呈現全部的資料
SELECT *, COUNT(*) 
FROM State_Info JOIN Farm_Info
ON State_Info.State_ID = Farm_Info.State_ID
GROUP BY State_Info.State_name,Farm_Info.Irrigation;

SELECT State_Info.State_namefg,COUNT(*) 
FROM State_Info JOIN Farm_Info
ON State_Info.State_ID = Farm_Info.State_ID
GROUP BY State_Info.State_name,Farm_Info.Irrigation;
