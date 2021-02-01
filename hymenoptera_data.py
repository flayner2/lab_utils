import pickle
import os
from typing import Union
from dotenv import load_dotenv
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request


def load_env_vars(env_path: str, key: str) -> Union[str, None]:
    load_dotenv(dotenv_path=env_path, override=True)

    return os.getenv(key)


# Environment variables configuration
env_path = os.path.join(
    os.path.expanduser("~/Documents/programming/python/projects/hymenoptera_data"),
    ".env",
)

# The ID and range of a sample spreadsheet.
SPREADSHEET_ID = load_env_vars(env_path, "SPREADSHEET_ID")
SPREADSHEET_RANGE = "data!A1:W72"

# If modifying these scopes, delete the file token.pickle.
SCOPES = ["https://www.googleapis.com/auth/spreadsheets.readonly"]


def main():
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
    else:

        for row in values[1:]:
            print(row)


if __name__ == "__main__":
    main()