#!/usr/bin/env python

###############
### Modules ###
###############
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.service import Service
import time
import argparse

#######################
### Argument parser ###
#######################
ap = argparse.ArgumentParser()
ap.add_argument("-s", "--smiles", required=True, help="smiles")
ap.add_argument("-o", "--outdir", required=True, help="outdir")
ap.add_argument("-p", "--phase", required=True, help="metabolism phase")
ap.add_argument("-g", "--gecko", required=True, help="geckodriver path")
args = vars(ap.parse_args())

##############
### Script ###
##############

options = Options()
options.add_argument('--headless')
options.set_preference("browser.download.folderList", 2)
options.set_preference("browser.download.manager.showWhenStarting", False)
options.set_preference("browser.download.dir", args['outdir'])
options.set_preference("browser.helperApps.neverAsk.saveToDisk", "text/csv")

service = Service(executable_path=args['gecko'])

driver = webdriver.Firefox(service=service, options=options)

driver.get("https://nerdd.univie.ac.at/gloryx/")

checkbox = driver.find_element(By.ID, "inputTextOption")

if not checkbox.is_selected():
    checkbox.click()

search_box = driver.find_element(By.ID, "input")
search_box.send_keys(args['smiles'])

dropdown = Select(driver.find_element(By.ID, "metabolism_phase"))

if args['phase'] == "All" :
    dropdown.select_by_visible_text("Phase 1 and phase 2 metabolism")
elif args['phase'] == 1 :
    dropdown.select_by_visible_text("Phase 1 metabolism")
elif args['phase'] == 2 :
    dropdown.select_by_visible_text("Phase 2 metabolism")

submit_button = driver.find_element(By.CSS_SELECTOR, "button.btn.btn-lg.btn-primary")
submit_button.click()

time.sleep(70)

csv_option = WebDriverWait(driver, 10).until(
    EC.presence_of_element_located((By.XPATH, "/html/body/div/main/header/section/div/div/div/div[3]/div[2]/div/a[2]"))
)

driver.execute_script("arguments[0].scrollIntoView(true);", csv_option)
driver.execute_script("arguments[0].click();", csv_option)

time.sleep(5)

driver.quit()
