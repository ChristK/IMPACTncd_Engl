# Installing R on Windows

#### Authors : Dr. Chris Kypridemos, Dr. Anna Head, Adithi R. Upadhya


- The installation procedure outlined here pertains to versions R 4.3.2. 

- The installation process for R has exhibited notable consistency throughout the years.

### Steps to install R on Windows

1. To install R navigate to [CRAN](https://cran.rstudio.com/)

2. Click **Download R for Windows** 

![](img/R_1.jpeg)


3. Once you click it, click on **install R for the first time**

![](img/R_2.jpeg)


4. Now click on **Download R-4.3.2 for Windows**

![](img/R_3.jpeg)


5. Double click on the downloaded *.exe* file which is the setup file, Click **Yes** (accept defaults)

![](img/R_4.jpeg)


6. Select the language needed, here we have chosen here **English**, the default one

![](img/R_5.jpeg)


7. Select default features, Click **Next** 

![](img/R_6.jpeg)


8. Click **Next** again to install using default settings and default location, you can also browse to another folder if required

![](img/R_7.jpeg)


9. Then click **Next**  to install the following components

![](img/R_8.jpeg)


10. Check **No**  for default start up options and click **Next** 

![](img/R_9.jpeg)


11. Click **Next** to select the start Menu folder (accept default)

![](img/R_10.jpeg)


12. Check on whichever options you need and click on **Next**

![](img/R_11.jpeg)


13. Wait for the installation process

![](img/R_12.jpeg)


14. Once completed this following window will show up, Click **Next**. Now we need to add paths of R in the system variables

15. Type Environment variables in search bar of your Windows system, click on **Edit Environment Variables**

![](img/Git_16.jpeg)


16. Click on New  and then add the variable name **PATH** and the value as *C:\\Program Files\\R-4.3.2\\bin* which is the path where bin folder or the Rscript.exe is stored e.g : *R4.3.2->bin->Rscript.exe*, click **Apply**

![](img/Git_17.jpeg)


19. To check Rscript works on command prompt or console or terminal, go to the search bar in Windows system, type **cmd** and then click on command prompt / console / terminal. Once it opens, type `Rscript` in terminal and you should be able to see this

![](img/Git_19.jpeg)

