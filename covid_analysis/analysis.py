import geopandas as gpd
from shapely.geometry import Point, Polygon
from shapely.affinity import translate, scale
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import re #regular expression

def plot_cases_on_world_map(dfs, date):
    def get_confirmed_cases(dfs, country): 
        confirmed = []
        dates = ["1-%s"%(i) for i in range(22,32)] + ["2-%s"%(i) for i in range(1,30)] + ["3-%s"%(i) for i in range(1,22)]
        for date in dates:
            if country == "United States of America":
                dfs[date]=dfs[date].replace(to_replace="US", value="United States of America", regex=False)
            val = dfs[date].loc[(dfs[date]['Country/Region'] == country)]['Confirmed'].sum()
            confirmed.append(val)
        confirmed = np.array(confirmed)
        return confirmed, np.array(dates)
    world = make_world_map()
    confirmed_cases = []
    for country in (world['name']):
        confirmed, dates = get_confirmed_cases(dfs, country)
        dfs[country] = pd.DataFrame()
        dfs[country]['dates'] = dates
        dfs[country]['COVID_cases']=confirmed
        confirmed_cases.append(dfs[country].loc[dfs[country]['dates']==date]['COVID_cases'].to_numpy())
    world['COVID_cases'] = np.array(confirmed_cases)
    world['COVID_density'] = (np.array([val[0] for val in confirmed_cases]))/world['pop_est']*1e6
    fig, ax = plt.subplots(nrows=2)
    ax1 = world.plot(column = 'COVID_cases', ax = ax[0], cmap = 'plasma', norm = matplotlib.colors.LogNorm(vmin = 0.9, vmax = 200000), legend = True, legend_kwds={'label': "COVID-19 cases", 'orientation': "vertical"})
    ax2 = world.plot(column = 'COVID_density', ax = ax[1], cmap = 'plasma', norm = matplotlib.colors.LogNorm(vmin = 0.9, vmax = 1000), legend = True, legend_kwds={'label':'Cases per million people', 'orientation': "vertical"})
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax1.set_title("COVID-19 confirmed cases %s"%(date))
    ax2.set_title("COVID-19 confirmed cases per million people %s"%(date))
    fig.savefig('world_maps/%s_update.png'%(date))
    #plt.show()
    plt.close()

def plot_cases_on_map(dfs, date):
    def get_confirmed_cases(dfs, state, state_full):
        confirmed = []
        dates = ['3-%s'%(i) for i in range(1,22)]
        for date in dates:
            dfs[date]=dfs[date].replace(to_replace='.*, %s'%(state), value='%s'%(state_full), regex=True)
            val = dfs[date].loc[(dfs[date]['Province/State'] == state_full)]['Confirmed'].sum()
            confirmed.append(dfs[date].loc[(dfs[date]['Province/State'] == state_full)]['Confirmed'].sum())
        confirmed = np.array(confirmed)
        return confirmed, np.array(dates)
    us = make_state_map()
    us = us.sort_values(by='Name')
    pop = pd.read_excel('us_pop.xlsx') #read in dataframe for state populations
    us = us.assign(pop_est=pop['pop_est'].to_numpy())
    us['geometry'][39] = translate(us['geometry'][39], 50, 5) #translate Hawaii
    us['geometry'][27] = scale(us['geometry'][27], 0.33, 0.33) #scale Alaska
    us['geometry'][27] = translate(us['geometry'][27], -70, -35) #translate Alaska
    confirmed_cases = []
    for abbr, state in zip(us['StateAbbr'], us['Name']):
        confirmed, dates = get_confirmed_cases(dfs, abbr, state)
        dfs[abbr] = pd.DataFrame()
        dfs[abbr]['dates'] = dates
        dfs[abbr]['COVID_cases']=confirmed
        confirmed_cases.append(dfs[abbr].loc[dfs[abbr]['dates']==date]['COVID_cases'].to_numpy())
    us['COVID_cases'] = np.array(confirmed_cases)
    us['COVID_density'] = (np.array([val[0] for val in confirmed_cases]))/us['pop_est']*1e6
    fig, ax = plt.subplots(nrows=2)
    ax1 = us.plot(column = 'COVID_cases', ax = ax[0], cmap = 'plasma', norm = matplotlib.colors.LogNorm(vmin = 0.9, vmax = 10000), legend = True, legend_kwds={'label': "Confirmed cases", 'orientation': "vertical"})
    ax1.set_xlim(-130, -64)
    ax1.set_ylim(22, 52)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_title("COVID-19 cases %s"%(date))
    ax2 = us.plot(column = 'COVID_density', ax = ax[1], cmap = 'plasma', norm = matplotlib.colors.LogNorm(vmin = 0.9, vmax = 1000), legend = True, legend_kwds={'label': "Cases per million people", 'orientation': "vertical"})
    ax2.set_xlim(-130, -64)
    ax2.set_ylim(22, 52)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_title("COVID-19 cases per million people %s"%(date))
    fig.savefig('state_maps/%s_updated.png'%(date))
    #plt.show()
    plt.close()
    
        
def make_state_map():
    usa = gpd.read_file('/Users/vahsengrouplaptop/Downloads/cb_2018_us_state_500k/cb_2018_us_state_500k.shp')
    usa = usa.drop(columns = ['STATEFP', 'STATENS', 'AFFGEOID', 'GEOID', 'LSAD'])
    usa.columns = ['StateAbbr', 'Name', 'Land_Area', 'Water_Area', 'geometry']
    usa = usa.loc[(usa['Name'] != 'Guam')]
    usa = usa.loc[(usa['Name']!='Commonwealth of the Northern Mariana Islands')]
    usa = usa.loc[usa['Name']!='American Samoa']
    usa = usa.loc[(usa['Name'] != 'United States Virgin Islands')]
    usa = usa.loc[usa['Name'] != 'District of Columbia']
    usa.index = [i for i in range(0,len(usa))]
    return usa

def make_world_map():
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    world=world.drop(columns = ['iso_a3'])
    return world

def plot_state(dfs, state, state_full, log):
    dat = [key for key in dfs.keys()]
    day = get_day(dat[len(dat)-1])
    dates = ['3-%s'%(i) for i in range(1,day+1)]
    dates_plt = []
    confirmed = []
    for date in dates:
        dfs[date]=dfs[date].replace(to_replace='.*, %s'%(state), value='%s'%(state_full), regex=True)
        val = dfs[date].loc[(dfs[date]['Province/State'] == state_full)]['Confirmed'].sum()
        if val > 0:
            dates_plt.append(date)
            confirmed.append(dfs[date].loc[(dfs[date]['Province/State'] == state_full)]['Confirmed'].sum())
    confirmed = np.array(confirmed)
    df = pd.DataFrame()
    df['confirmed'] = confirmed
    x = np.array([i for i in range(0,len(confirmed))])
    if len(df.loc[df['confirmed']>0]) > 1:
        fit = np.polyfit(x, np.log(confirmed), 1)
    else:
        fit = np.array([0,0])
        print("Fit didn't converge: fit set to 0")
    a = 1/fit[0]
    b = -fit[1]/fit[0]
    x = np.array([i for i in range(0,len(confirmed)+2)])
    
    for i in range(0, len(dates)):
        if dates[i] == dates_plt[0]:
            shift = i

    plt.plot(['3-%s'%(i) for i in range(1,32)], [0 for i in range(0,31)], 'o', markersize = 0)
    plt.plot(dates_plt, confirmed, 'o', color = 'tab:blue')
    plt.plot(x+shift, np.exp(fit[0]*x+fit[1]), '--', color = 'red', label = 'exponential fit, doubling time = %s days'%(float('%.3g' % a)))
    locs, labels = plt.xticks()
    plt.xticks(np.arange(0, 31, step=3))
    plt.ylabel("COVID-19 Cases")
    plt.title("Confirmed Cases %s"%(state_full))
    if log == "log":
        plt.yscale("log")
    plt.legend()
    #plt.show()
    plt.savefig('state_fits/%s_3-%s.png'%(state_full, day))
    plt.close()
    return confirmed.max(), fit, len(confirmed)

def plot_fit(dfs, country, log, color): #choose lin or log
    date = [key for key in dfs.keys()]
    day = get_day(date[len(date)-1])
    fit= fit_exponential_to_march_data(dfs, country)
    a = 1/fit[0]
    b = -fit[1]/fit[0]
    confirmed = []
    dates = ['3-%s'%(i) for i in range(1,day+1)]
    for date in dates:
        if country == "World":
            confirmed.append(dfs[date]['Confirmed'].sum())
        else:
            confirmed.append(dfs[date].loc[dfs[date]['Country/Region'] == country]['Confirmed'].sum())
    for i in range(0,len(dates)):
        if i == 0:
            plt.plot(dates[i], confirmed[i], 'o', color = '%s'%(color),  markersize = 3, label = '%s Confirmed Cases. Fit doubling time = %s days'%(country, float('%.3g' % a)))
        else:
            plt.plot(dates[i], confirmed[i], 'o', color = '%s'%(color),  markersize = 3)
    if log == "log":
        x = np.array([i for i in range(0,31)])
        plt.plot(x, np.exp((x-b)/a), '--', color = '%s'%(color))
        plt.plot(dates + ['3-%s'%(i) for i in range(day+1, 32)], [0 for i in range(0,31)], 'o', markersize = 0)
        plt.yscale("log")
        #plt.ylim(1, 500000)
    else:
        x = np.array([i for i in range(0,day+2)])
        plt.plot(x, np.exp((x-b)/a), '--', color = '%s'%(color))
        plt.plot(dates + ['3-%s'%(i) for i in range(day+1, day+3)], [0 for i in range(0,day+2)], 'o', markersize = 0)
        #plt.ylim(0,1500)
    locs, labels = plt.xticks()
    plt.xticks(np.arange(0, len(x), step=3))
    #plt.xlabel("date")
    plt.ylabel("COVID-19 Cases")
    plt.legend()
    #plt.savefig('country_fits/%s_3-%s.png'%(country, day))
    #plt.close()
    
    
def plot_confirmed_cases(dfs):
    date = [key for key in dfs.keys()]
    day = get_day(date[len(date)-1])
    confirmed_US = []
    confirmed_italy = []
    confirmed = []
    confirmed_china = []
    for key in dfs.keys():
        confirmed_US.append(dfs[key].loc[dfs[key]['Country/Region'] == "US"]['Confirmed'].sum())
        confirmed_italy.append(dfs[key].loc[dfs[key]['Country/Region'] == "Italy"]['Confirmed'].sum())
        confirmed_china.append(dfs[key].loc[dfs[key]['Country/Region'] == "China"]['Confirmed'].sum())
        confirmed.append(dfs[key]['Confirmed'].sum())
    plt.plot([key for key in dfs.keys()], confirmed, 'o', color = 'black', markersize = 3, label = 'Confirmed Cases Worldwide')
    plt.plot([key for key in dfs.keys()], confirmed_italy, 'o', color = 'red', markersize = 3, label = 'Confirmed Cases Italy')
    plt.plot([key for key in dfs.keys()], confirmed_china, 'o', color = 'green', markersize = 3, label = 'Confirmed Cases China')
    plt.plot([key for key in dfs.keys()], confirmed_US, 'o', color = 'blue', markersize = 3, label = 'Confirmed Cases US')
    fit= fit_exponential_to_march_data(dfs, "US")
    x = np.array([i for i in range(1,32)])
    plt.plot(range(40, 40+len(x)), np.exp(fit[0]*x + fit[1]), '--', color = 'tab:blue', label = "US projected number of cases assuming exponential growth")
    plt.plot([key for key in dfs.keys()] + ["3-%s"%(i) for i in range(day+1,32)], [0 for i in range(0,len([key for key in dfs.keys()] + ["3-%s"%(i) for i in range(day+1,32)]))], 'o', markersize = 0) #Expand x range to include future dates
    #plt.plot([len(dfs.keys())-6,len(dfs.keys())+7], [confirmed[len(confirmed)-6],confirmed[len(confirmed)-6]], '--', color = 'black', linewidth = 0.5)
    #plt.plot([len(dfs.keys())-6,len(dfs.keys())], [confirmed_italy[len(confirmed_italy)-6],confirmed_italy[len(confirmed_italy)-6]], '--', color = 'red', linewidth = 0.5)
    #plt.plot([len(dfs.keys())-6,len(dfs.keys())+4], [confirmed_china[len(confirmed_china)-6],confirmed_china[len(confirmed_china)-6]], '--', color = 'green', linewidth = 0.5)
    locs, labels = plt.xticks()
    plt.xticks(np.arange(0, len([key for key in dfs.keys()] + ["3-%s"%(i) for i in range(day+1,32)]), step=3))
    plt.title("Confirmed COVID-19 Cases over time")
    plt.legend()
    plt.ylim(0,350000)
    #plt.xlabel("month-day")

def projection_from_death_rate(dfs, state, death_rate):
    doubling_rate = 3.5
    date = [key for key in dfs][len(dfs)-1]
    date1 = [key for key in dfs][len(dfs)-2]
    deaths = dfs[date].loc[dfs[date]['Province/State'] == "%s"%(state)]['Deaths'].sum()
    dt = 17 #days between when person got the virus and when they died
    date0 = [key for key in dfs][len(dfs)-1-dt]
    cases = dfs[date0].loc[dfs[date0]['Province/State'] == "%s"%(state)]['Confirmed'].sum()
    
    
def fit_exponential_to_march_data(dfs, country):
    date = [key for key in dfs.keys()]
    day = get_day(date[len(date)-1])
    confirmed = []
    #dates = ['3-%s'%(i) for i in range(1,17)]
    dates = ['3-%s'%(i) for i in range(1,day+1)]
    for date in dates:
        if country == "World":
            confirmed.append(dfs[date]['Confirmed'].sum())
        else:
            confirmed.append(dfs[date].loc[dfs[date]['Country/Region'] == "%s"%(country)]['Confirmed'].sum())
    #x = np.linspace(0,15,16)
    x = np.linspace(0,day-1,day)
    confirmed = np.array(confirmed)
    fit = np.polyfit(x, np.log(confirmed),1)
    return fit
    
def make_dataframes(date):
    day = get_day(date)
    dates = ["01-%s-2020"%(i) for i in range(22,32)] + ["02-0%s-2020"%(i) for i in range(1,10)] + ["02-%s-2020"%(i) for i in range(10,30)] + ["03-0%s-2020"%(i) for i in range(1,10)] + ["03-%s-2020"%(i) for i in range(10,day+1)]
    dfs = {}
    for date in dates:
        dfs[date] = pd.read_csv("csse_covid_19_daily_reports/%s.csv"%(date))
        dfs[date]=dfs[date].replace(to_replace="Cote d'Ivoire", value="Côte d'Ivoire", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Mainland China", value="China", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Taiwan*", value="Taiwan", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Congo (Kinshasa)", value="Dem. Rep. Congo", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Congo (Brazzaville)", value="Congo", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Mayotte", value="France", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Holy See", value="Italy", regex=False)
        dfs[date]=dfs[date].replace(to_replace="The Bahamas", value="Bahamas", regex=False)
        dfs[date]=dfs[date].replace(to_replace="The Gambia", value="Gambia", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Guernsey", value="United Kingdom", regex=False)
        dfs[date]=dfs[date].replace(to_replace="occupied Palestinian territroy", value="Israel", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Jersey", value="United Kingdom", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Republic of the Congo", value="Congo", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Korea, South", value="South Korea", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Equatorial Guinea", value="Eq. Guinea", regex=False)
        dfs[date]=dfs[date].replace(to_replace="San Marino", value="Italy", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Dominican Republic", value="Dominican Rep.", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Central African Republic", value="Central African Rep.", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Andorra", value="France", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Malta", value="Italy", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Bosnia and Herzegovina", value="Bosnia and Herz.", regex=False)
        dfs[date]=dfs[date].replace(to_replace="North Macedonia", value="Macedonia", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Singapore", value="Malaysia", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Monaco", value="France", regex=False)
        dfs[date]=dfs[date].replace(to_replace="Bahrain", value="Qatar", regex=False)
        dfs[date] = dfs[date].rename(columns={"Country_Region": "Country/Region", "Province_State": "Province/State"})
    dfs = { key.replace('-2020', ''):value for key, value in dfs.items() } #remove the -2020's from the key
    for i in range(0,10):
        dfs = { key.replace('0%s'%(i), '%s'%(i)):value for key, value in dfs.items() } #removes 0's in the date keys
    return dfs

def generate_report(dfs):
    plt.subplot(3,1,1)

    plot_confirmed_cases(dfs)

    plt.subplot(3,1,2)
    #plot_fit(dfs, "World", "lin", "black")
    plot_fit(dfs, "US", "lin", "blue")
    #plot_fit(dfs, "Germany", "lin", "green")
    #plot_fit(dfs, "France", "lin", "gold")
    #plot_fit(dfs, "Italy", "lin", "red")

    plt.subplot(3,1,3)
    #plot_fit(dfs, "World", "log", "black")
    plot_fit(dfs, "US", "log", "blue")
    #plot_fit(dfs, "Germany", "log", "green")
    #plot_fit(dfs, "France", "log", "gold")
    #plot_fit(dfs, "Italy", "log", "red")

    plt.show()

def get_day(date):
    day = re.findall(r'\d+', date) # takes out the integers from the date string
    day = int(day[1])
    return day

dfs = make_dataframes('3-25')
generate_report(dfs)



'''
def plot_cases_on_map_projection(log):
    us = make_state_map()
    cases = []
    fits = []
    num_days = []
    for abbr, state in zip(us['StateAbbr'], us['Name']):
        confirmed, fit, num_day = plot_state(abbr, state, "linear")
        cases.append(confirmed)
        fits.append(np.exp(fit[0]*(num_day)+fit[1]))
        num_days.append(num_day)
    fig, ax = plt.subplots(nrows=2)
    us['num_days'] = np.array(num_days)
    if log == "log":
        us['COVID_cases'] = np.log10(np.array(cases)+1)
        us['Proj_cases'] = np.log10(np.array(fits)+1)
        ax1 = us.plot(column = 'COVID_cases', ax = ax[0], cmap = 'plasma', vmin = 0, vmax = 4, legend = True, legend_kwds={'label': "log(COVID-19 cases)", 'orientation': "horizontal"})
        ax2 = us.plot(column = 'Proj_cases', ax = ax[1], cmap = 'plasma', vmin = 0, vmax = 4, legend = True, legend_kwds={'label': "log(Projected cases tomorrow)", 'orientation': "horizontal"})
    else:
        us['COVID_cases'] = np.array(cases)
        us['Proj_cases'] = np.array(fits)
        ax1 = us.plot(column = 'COVID_cases', ax = ax[0], cmap = 'plasma', legend = True, legend_kwds={'label': "COVID-19 cases", 'orientation': "horizontal"})
        ax2 = us.plot(column = 'Proj_cases', ax = ax[1], cmap = 'plasma', legend = True, legend_kwds={'label': "Projected cases tomorrow", 'orientation':"horizontal"})
    ax1.set_xlim(-130, -64)
    ax1.set_ylim(24, 52)
    ax2.set_xlim(-130, -64)
    ax2.set_ylim(24, 52)
    fig.show()
'''
