

function click_collapse_samstat() {
   var samstat_block = document.getElementsByClassName("samstat_block"); 
    if (samstat_block.style.display == "block") {
        samstat_block.style.display = "none";
    } else {
        samstat_block.style.display = "block";
    }
}
