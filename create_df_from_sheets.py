import pickle
import os
import pandas as pd
from typing import Union
from dotenv import load_dotenv
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request


def load_env_vars(key: str) -> Union[str, None]:
    load_dotenv(override=True)

    return os.getenv(key)


# The ID and range of a sample spreadsheet.
SPREADSHEET_ID = load_env_vars("SPREADSHEET_ID")
SPREADSHEET_RANGE = load_env_vars("SPREADSHEET_RANGE")

# If modifying these scopes, delete the file token.pickle.
SCOPES = ["https://www.googleapis.com/auth/spreadsheets.readonly"]


def load_data() -> Union[pd.DataFrame, None]:
    """Shows basic usage of the Sheets API.
    Prints values from a sample spreadsheet.
    """
    creds = None
    # The file token.pickle stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    if os.path.exists("token.pickle"):
        with open("token.pickle", "rb") as token:
            creds = pickle.load(token)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file("credentials.json", SCOPES)
            creds = flow.run_local_server(port=0)
        # Save the credentials for the next run
        with open("token.pickle", "wb") as token:
            pickle.dump(creds, token)

    service = build("sheets", "v4", credentials=creds)

    # Call the Sheets API
    sheet = service.spreadsheets()
    result = (
        sheet.values()
        .get(spreadsheetId=SPREADSHEET_ID, range=SPREADSHEET_RANGE)
        .execute()
    )
    values = result.get("values", [])

    if not values:
        print("No data found.")

        return None
    else:
        # Change this to fit the shape of your own data
        data = pd.DataFrame(data=values[1:], columns=values[0])

        return data


if __name__ == "__main__":
    load_data()
