#!/bin/bash

read -p "Type Protocol name: " arg_protocol_name
echo "$arg_protocol_name"
read -p "Type Arm name: " arg_arm_name
echo "$arg_arm_name"
read -p "Type filename to save data in: " arg_filename
echo "$arg_filename"
echo "Type mysql db password: "
read -s db_pass

# query="select b.id_genealogy, q.volume, s.date, q.notes from quantitative_measure as q, bio_mice as b, measurements_series as s where q.id_mouse in (select id_mouse from mice_has_arms where id_protocols_has_arms in (select id from protocols_has_arms where id_protocol in (select id from protocols where name=\"$arg_protocol_name\") and id_arm in (select id from arms where name like \"%$arg_arm_name%\"))) and q.id_series in (select id_series from measurements_series where id_type=2) and q.id_mouse=b.id and q.id_series=s.id_series order by b.id_genealogy, s.date asc;"
query="select b.id_genealogy, g.name, p.name, a.name, q.volume, s.date, mha.start_date, mha.expected_end_date, mha.end_date, q.notes from
 quantitative_measure as q, bio_mice as b, measurements_series as s, groups as g, mice_has_arms as mha, protocols as p, arms as a
  where q.id_mouse in (select id_mouse from mice_has_arms where id_protocols_has_arms in (select id from protocols_has_arms where id_protocol in (select id from protocols where name=\"$arg_protocol_name\") 
  and id_arm in (select id from arms where name like \"%$arg_arm_name%\")))
  
  and q.id_series in (select id_series from measurements_series where id_type=2) 
  and q.id_mouse=b.id and q.id_series=s.id_series and b.id_group=g.id and mha.id_mouse=b.id 
  and p.id in (select id from protocols where name = \"$arg_protocol_name\") 
  and a.id in (select id from arms where name like \"%$arg_arm_name%\") order by b.id_genealogy, s.date asc"
# echo "$query"
echo -e "Measures of all biomice treated with ${arg_protocol_name} for arm ${arg_arm_name} \n" >> ./"${arg_filename}"

mysql -u root -p$db_pass -D caxeno -e "$query;" >> ./"${arg_filename}"

echo "Output saved in $arg_filename"
