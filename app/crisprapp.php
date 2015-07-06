<?php

// Define the variables
$start = $_POST['start'];
$end = $_POST['end'];



// Main process
if ($start && $end):

   $tmp = "/tmp/result.out";
   
   // read the file:
   $handle = fopen("./oligos/results.txt", "r");
   
   if ($handle) {
      while (($line = fgets($handle)) !== false) {
        // process the line read.
        $split = preg_split("/[\t,]/", $line); 
    }

    fclose($handle);
} else {
    // error opening the file.
} 

?>

<h3><font color="#00B85C"><b><u>CRISPR</u><sub>App</sub> <u>Results</b></u></font></h3>

<?php

include "/tmp/result.out";


else: ?>

<font color="red"><b>Sequence missing!</b></font> 

<?php endif; 


