/**
 *Batch macro to merge data from Operetta
 * experiments
 *
 *
 */

    //break out here
    //setBatchMode("exit and display");
    //stop

//Check that all prepared
{
    //waitForUser("Files must be arranged as T1->T2->T3(z)->Tub sequence. Proceed?");
}

//INPUT VARIABLES
{
  // Run options:
  test = false;
  test_size = 200; //cycles to test
  start = 1; //file to start from

  //Type of projection
  zproj_array = newArray("Max Intensity","Average Intensity");
  zproj_ind_array = newArray("max","ave");
  project_type = 0;

  //imaging settings
  channels = 4;
  zstacks = 3;
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
    resDir = getDirectory("Save data to:");
}

// Get file list
list = getFileList(dataDir);

//Processing loop
step = channels*zstacks;
firstfile = (start-1)*step;
if (test) {
  files2process = firstfile+test_size;
}
else{
  files2process = list.length;
}

for (fov=firstfile; fov<files2process; fov+=step)
{
  //Get name of image i
  {
      fileNameFull = list[fov];
      //dotIndex = lastIndexOf(fileNameFull, ".");
      cutIndex = lastIndexOf(fileNameFull, "-");
      if (cutIndex!=-1);
      fileName = substring(fileNameFull, 0, cutIndex); // name of file w/o extension
  }

  //print status
  current_file = (fov+step)/step;
  files_total =  list.length/step;
  complete = d2s((current_file/files_total)*100,0);
  print("---");
  print("File " + d2s(current_file,0) + " of " + files_total);
  print(fileName);
  print("["+complete+"%]");

  //Open all images
  step2 = zstacks;
  for (channel=0; channel<step; channel+=step2) {
      //open z slices for each channel
      for (z=0; z<zstacks; z+=1) {
        path2slice = dataDir+list[fov+channel + z];
        open(path2slice);
      }
  }

  //Move images to Stack
  {
    run("Images to Stack", "name=Stack title=[] use");
    //Convert to hyperstack
    run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices="+zstacks+" frames=1 display=Grayscale");
    //Make Z projection
    run("Z Project...", "projection=["+zproj_array[project_type]+"]");
    rename("Result");
  }

  //Save individual Channels
  for (channel=1; channel<channels+1; channel+=1) {
    selectWindow("Result");
    run("Duplicate...", "title=ch"+channel+" duplicate channels="+channel);
    save(resDir+fileName+"_ch"+channel+"_"+zproj_ind_array[project_type]+".tif");
  }

  //Clean
  {
      //Closes all images:
      while (nImages>0) {
          selectImage(nImages);
          close();
      }
  }
}
//print status
print("==================================");
print("===DONE===========================");
