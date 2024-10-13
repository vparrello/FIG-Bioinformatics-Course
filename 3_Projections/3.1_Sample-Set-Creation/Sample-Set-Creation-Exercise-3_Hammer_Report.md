# Sample Set Creation: Generating Hammer Reports

## Objective
Learn how to generate hammer reports for each of the samples in your dataset, which will be crucial for further analysis and machine learning model creation.

## Materials
[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)
[NCBI](https://www.ncbi.nlm.nih.gov/)
[SRAToolkit Download FASTA Documentation](https://www.ncbi.nlm.nih.gov/books/NBK242621/)

FIG-Bioinformatics-Course/
├── 3_Projections
    └── 3.1_Sample-Set-Creation/
        └── Sample-Set-Creation Exercise-3_Hammer_Report.md (you are here)


## Exercise

1. Review the concept of hammer reports. Ask Grimoire to explain what a hammer report is in the context of metagenomics and why it's important for our machine learning project.

2. Ensure that you have the FIGfam database installed and accessible. If not, follow the installation instructions provided in the course materials or ask for assistance.

3. Familiarize yourself with the FIGARO tool, which we'll use to generate the hammer reports. Ask Grimoire to explain how FIGARO works and its key features.

4. For each sample in your dataset (both control and diseased groups), you'll need to generate a hammer report. Here's a general outline of the process:

   a. Navigate to the directory containing your sample data.
   b. Use the FIGARO tool to generate a hammer report for each sample. The basic command structure is:

      ```
      figaro -i <input_file> -o <output_file> -d <figfam_database>
      ```

      Replace `<input_file>` with your sample file name, `<output_file>` with the desired name for your hammer report, and `<figfam_database>` with the path to your FIGfam database.

5. Write a simple shell script that automates this process for all your samples. Ask Grimoire for help if you're not familiar with shell scripting.

6. Run your script and ensure that you have a hammer report for each of your samples.

7. Examine a few of the hammer reports to familiarize yourself with their structure and content. Ask Grimoire to help you interpret key sections of the report.

8. Organize your hammer reports in a clear directory structure, separating reports for control and diseased samples.

9. Reflect on how these hammer reports will be used in the next steps of your machine learning project. What kind of information do they provide that will be useful for classification?

10. (Optional) If time allows, ask Grimoire about different ways to visualize or summarize the data from your hammer reports. This could provide interesting insights into your dataset before you move on to the machine learning phase.

Remember, generating hammer reports for a large number of samples can be time-consuming. Plan accordingly and consider running your script overnight if you have a particularly large dataset.



Own Lesson: Get a hammer report from each of the samples