<?php
$hostname = $_POST["hostname"];
$username = $_POST["username"];
$password = $_POST["password"];
$myquery  = $_POST["myquery"];
$link           = mysql_connect("$hostname","$username","$password");
if(!$link) {print("i'm dead!\n" . mysql_error() . "\n");};
$explode_result = explode("_pck_",$myquery);
$L              = count($explode_result);
for ($ii=0;$ii<$L;$ii=$ii+1)
{
  $query = $explode_result[$ii];
  $query = str_replace("_apice_","'",$query);
  $query = str_replace("_spazio_"," ",$query);
  $query = str_replace("_slash_","\\",$query);
  $query = str_replace("_virgo_","\"",$query);
  $query = str_replace("_canc_","#",$query);
  $query = str_replace("_perc_","%",$query);
  $query = str_replace("_and_","&",$query);
  $query = str_replace("_at_","@",$query);
  $query = str_replace("_equals_","=",$query);
  $query = str_replace("_period_",".",$query);
  $query = str_replace("_star_","*",$query);
  $query = str_replace("_colon_",":",$query);
  $query = str_replace("_comma_",",",$query);



  $risultato = mysql_query($query);
  if(is_resource($risultato)) {
   while ($riga = mysql_fetch_row($risultato))
     {
      $N = count($riga,1);
      for ($jj=0;$jj<$N;$jj=$jj+1)
      {
       print("$riga[$jj]\n");
      }
      print("\n");
     }
   }
   print("\r");
}

?>
