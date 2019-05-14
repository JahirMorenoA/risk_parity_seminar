# Risk Parity Seminar

To create a risk-parity-seminar image run

```
cd $(pwd)
docker build -t risk-parity-seminar .
docker run --rm -e PASSWORD=book -p 8787:8787 -v $(pwd)/seminar:/home/rstudio/seminar risk-parity-seminar
```

Then go to `localhost:8787` in the browser. 

Within Rstudio use rstudio/seminar when asked for user/password and go to /home/rstudio/seminar/



