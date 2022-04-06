/**
 *Batch macro to make low res tile of whole well of 20x data (1 img per well)
 * from Operetta datasets
 */

    //break out here
    //setBatchMode("exit and display");
    //stop

// Run options:
{
  test = false;
  //Aquisition parameters
  channel_number = 4;

  //FOVs per well
  //img_N = 25;

  //Montage options
  mont_opt = "columns=5 rows=5 scale=1 last=25 border=1"

  //crop width
  w = 200 //pxl
}

//Prepare workspace:
{
    run("Close All");
    print("\\Clear");
    run("Clear Results");
}

setBatchMode(true);

// Get paths
{
    dataDir = getDirectory("Process next files:");
    //print(dataDir)
    //dataDir = 'D:/Temp/APID/APID_5_Screen2/Morpho/8_SignatureAnalysis/CellImages/Output/Raw/';
    resDir = getDirectory("Save data to:");
    //resDir = 'D:/Temp/APID/APID_5_Screen2/Morpho/8_SignatureAnalysis/CellImages/Output/Raw/';
}

// Get file list
list = getFileList(dataDir);

//get cellN
cellN = lengthOf(list)/channel_number

//Get name of the donor
{
  path = dataDir;
  folder = File.getName(path);
}

for (channel=1; channel<=channel_number; channel+=1)
{
  //pos of firt file
  start = channel-1;

  //Prepare empty stack
  newImage("Stack"+channel, "16-bit black", w, w, cellN);

  //load all files for the channel into Stack
  for (cell = 0; cell < cellN; cell+=1)
  {
    //Path to img i
    img_path = dataDir+list[cell*channel_number + channel-1];
    //open image & move it to the stack
    open(img_path);
    rename("img");
    run("Select All");
    run("Copy");
    selectWindow("Stack"+channel);
    Stack.setSlice(cell+1);
    run("Paste");
    selectWindow('img');
    run("Close");
  }

  //Make montage
  {
    selectWindow("Stack"+channel);
    run("Make Montage...", mont_opt);
    rename("Montage"+channel);
    selectWindow("Stack"+channel);
    run("Close");
  }

}

//Merge Channels
/** Colors:
  c1 = red
  c2 = green
  c3 = blue
  c4 = gray
  c5 = cyan
  c6 = magenta
  c7 = yellow
*/
run("Merge Channels...", "c1=Montage3 c2=Montage2 c5=Montage1 c6=Montage4 create");
rename("Result");

//Save Results
save(resDir+folder+'_montage.tif');

//Clean
{
    //Closes all images:
    while (nImages>0) {
        selectImage(nImages);
        close();
    }
}
//print status
print("==================================");
print('===DONE, donor = '+folder+'===========================');
