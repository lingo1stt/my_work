show databases;
use group_datascience;
SELECT * from Badges 
limit 100;
describe Badges ;

SELECT * from Comments  
limit 100;
describe Comments  ;

SELECT * from  PostHistory 
limit 100;
describe PostHistory  ;

SELECT * from  PostLinks  
limit 100;
describe PostLinks ;

SELECT * from  Posts 
limit 50;
describe Posts ;

SELECT * from  Tags  
limit 50;
describe Tags  ;

SELECT * from  Users 
limit 50;
describe Users  ;
##
SELECT Id ,reputation,CreationDate ,LastAccessDate ,Location ,Views ,UpVotes ,DownVotes  from Users

SELECT * from  Votes 
limit 50;
describe Votes  ;


##############################
#use users 
#############################
SELECT * from  Users 
limit 50;
describe Users  ;
##檢查row 數量
SELECT count(*) from Users ;
##reputation column
SELECT MAX(reputation), MIN(reputation) from Users ; 
SELECT AVG(reputation) from Users ; 

###取前100/後100出來分析
with percent as
(SELECT Reputation , Views,UpVotes, DownVotes,DisplayName  from Users 
ORDER BY  Reputation DESC
limit 100)
select AVG(Views), AVG(UpVotes), AVG(DownVotes) from percent;

with percent as
(SELECT Reputation , Views,UpVotes, DownVotes,DisplayName  from Users 
ORDER BY  Reputation ASC 
limit 100)
select AVG(Views), AVG(UpVotes), AVG(DownVotes) from percent;

SELECT * from Users 
WHERE Reputation  = 28918;
###date column

WITH TT AS (
  SELECT DATEDIFF(LastAccessDate, CreationDate) AS TIME_DIFF, DisplayName 
  FROM Users
)
SELECT DisplayName, TIME_DIFF 
FROM TT
WHERE TIME_DIFF = 3614;


### find the guy
SELECT * from Users 
where DisplayName  = 'rapaio' or 
DisplayName  = 'Sean Owen';


####location info
SELECT DisplayName, Location, from Users 

###about me
SELECT DisplayName, AboutMe  from Users 

C:\Users\lingo1st\OneDrive\桌面\groupwork_2.sql

### find the most popular post and its information 
#先把USER與POST TABLE進行連結
SELECT Users.Id,Posts.CreationDate ,Users.DisplayName ,Users.Location ,Posts.ViewCount ,Posts.Body ,Posts.Title ,Posts.Tags ,Posts.AnswerCount,Posts.CommentCount  FROM Users JOIN
Posts ON Users.Id =Posts.OwnerUserId ;





