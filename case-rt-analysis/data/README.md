# hMPXV-1 case data

## Preparing case data

Raw case counts by location are downloaded from Our World in Data's [monkeypox github repository](https://github.com/owid/monkeypox).

These counts are downloaded and filtered to locations and a time period of interest with the script `case-data-prep.py`:

```
python3 case-data-prep.py
```
